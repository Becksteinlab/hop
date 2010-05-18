# $Id$
# Hop --- a framework to analyze solvation dynamics from MD simulations
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Based on a definition of grid sites, convert a molecular
dynamics trajectory into a trajectory of site hops.

You will also need the following modules to create the input for HoppingTraj:
hop.sitemap, MDAnalysis
"""
import numpy
import MDAnalysis

import hop.constants
from hop.constants import SITELABEL
import hop.utilities
from hop.utilities import msg,set_verbosity
from hop import SelectionError

def totaltime(trajectory,targetunit='ps',sourceunit='AKMA'):
    """Returns the total trajectory time from the DCDReader object."""
    return trajectory.numframes * delta_t(trajectory,targetunit,sourceunit)

def delta_t(trajectory,targetunit='ps',sourceunit='AKMA'):
    """Returns the length of the time interval between two trajectory frames."""
    return trajectory.delta * trajectory.skip_timestep * hop.constants.get_conversion_factor('time',sourceunit,targetunit)

class HoppingTrajectory(object):
    """Provides a time-sequence of sites visited by individual molecules,
    called a 'hopping trajectory' because the molecules hop between
    sites. Their coordinates are mapped to site labels, which have been defined
    on a grid previously (using hop.sitemap).

    :Output format:

    For simplicity and code reusal this is again a dcd with the site as the
    x-coordinate; the y coordinate is set to the 'orbit site', i.e. it records
    the site the particle was last at for as long as it does not enter a new
    site. It describes the site in whose 'basin of attraction' the particle
    orbits. Note, however, that the transition to a new site is still counted
    as belonging to the previous site (which is arguably incorrect); the
    hop.graph module, however, does a proper analysis, which is cannot be done
    here for efficieny reasons. The z field is unused at the moment and set to
    0.

    :Attributes:

    ts               MDAnalysis.Timestep object
    numframes        number of frames in hopping trajectory
    group            AtomGroup of atoms that are tracked


    :Methods:

    ## [start:stop]       object can be used as an iterator over the
    ##                    hopping trajectory (disabled du to problems when doing random
    ##                    access on large dcds; either a bug in DCDReader or python)
    next()                advances time step in the hopping trajectory
    map_dcd()             iterator that updates the ts and maps the trajectory
                          coordinates to site labels
    _map_next_timestep()  map next coordinate trajectory step to hopping time step
    _read_next_timestep() read next timestep from hopping trajectory


    write()              write the hopping trajectory to a dcd file + psf
    write_psf()          write a dummy psf for visualization
    """

    def __init__(self,trajectory=None,group=None,density=None,
                 filename=None,hopdcd=None,hoppsf=None,fixtrajectory=None,verbosity=3):
        """Converts a trajectory into a hopping trajectory, using a sitemap as an index for sites.

        >>> h = HoppingTrajectory(trajectory=DCDReader,group=AtomGroup,density=Density,
                                  fixtrajectory=<dict>,verbosity=3)
        >>> h = HoppingTrajectory(filename=<name>)

        Create from a coordinate trajectory of a group of atoms and a site map:

          u = MDAnalysis.Universe(psf,dcd)
          water = u.selectAtoms('name OH2')
          h = HoppingTrajectory(trajectory=u.dcd,group=water,density=water_density)

        Load from a saved hopping trajectory (in dcd format with dummy psf)

          h = HoppingTrajectory(hopdcd='hops.dcd',hoppsf='hops.psf')

        :Arguments:

        trajectory       MDAnalysis.dcd trajectory instance
        group            MDAnalysis.group instance
        density          grid3Dc.Grid instance with sitemap set

        hopdcd           dcd written by write()
        hoppsf           psf written by write() (or write_psf())
        filename         or simply provide one filename prefix for psf and dcd

        fixtrajectory    dictionary with attributes of a dcd object and new
                         values; used to provide correct values after using
                         a catdcd-generated trajectory (hack!), e.g.
                         fixtrajectory = {'delta':10.22741474887299}
        
        verbosity        show status messages for >= 3
        """
        self.verbosity = verbosity
        set_verbosity(self.verbosity)

        if not (trajectory is None or group is None or density is None):
            self.traj  = trajectory        # MDAnalysis.Universe.trajectory
            self.tgroup = group            # atom selection for trajectory
            if not isinstance(self.tgroup,MDAnalysis.AtomGroup.AtomGroup):
                raise TypeError('group must be a <AtomGroup>, eg MDAnalyis.Universe.selectAtoms().')
            if isinstance(fixtrajectory,dict):
                for attr,val in fixtrajectory.items():
                    if not hasattr(trajectory,attr):
                        raise AttributeError('fixtrajectory: dcd object does not have attribute "'\
                                             +str(attr)+'"')
                    trajectory.__dict__[attr] = val            
            self.totaltime = totaltime(trajectory,'ps')
            self.traj.rewind()             # make sure to start from frame 0
            self._GD = density             # sitemap.Density object
            self.map   = self._GD.map                  # map of sites
            self.edges = self._GD.edges                # N+1 edges of bins 
            self.dedges = map(numpy.diff,self.edges)   # N bin widths
            try:
                if not self._GD.grid.shape == self.map.shape:
                    raise ValueError            
            except (AttributeError,ValueError):
                raise ValueError("The density object must have its site map computed.")
            Dmap = numpy.rank(self.map)
            coord = numpy.asarray(self.tgroup.coordinates())
            Natoms,D = coord.shape
            if not D == Dmap:
                raise ValueError("Coordinates and map have different dimensions.")
            # NOTE:
            # Any count outside the histogram becomes 'outlier' so
            # one should take care to choose a large enough map for the region
            # of interest. See _coord2hop().
            self.buffered_map = SITELABEL['outlier'] * \
                                numpy.ones(tuple(numpy.asarray(self.map.shape) + 2))

            # Here we commit to writing a DCD hopping trajectory:
            self.ts = MDAnalysis.DCD.Timestep(Natoms)   # empty time step for N atoms
            self.ts.frame = self.traj.ts.frame          # current frame
            numlabels = float(self.map.max() - self.map.min() + 2) # naive... but not crucial
            # fake unit cell for visualization
            # Layout of DCD unitcell is [A, alpha, B, beta, gamma, C] (sic!)
            self.ts._unitcell = numpy.array((numlabels,90, numlabels,90, 90,1),dtype=numpy.float32)        
            # current hopping trajectory frame is in ts._pos[]
            # _pos = numpy.empty(coord.shape)   # x=site label y=s(t)==0?s(t-1):s(t)  z=0
            self.numframes = self.traj.numframes    # total numer of frames
            self._init_coord2hop()                  # init for _map_next_timestep()
            self._map_next_timestep()               # initialize with first timestep
            self.hoptraj = None                     # no hopping trajectory available
        elif not (hopdcd is None or hoppsf is None) or filename is not None:
            # read from dcd
            try:
                self.traj,self.tgroup,self.map,self.edges,self.dedges
            except AttributeError:
                self.traj,self.tgroup,self.map,self.edges,self.dedges = [None] * 5
            if filename is not None:
                hoppsf = self.filename(filename,'psf')
                hopdcd = self.filename(filename,'dcd')
            u = MDAnalysis.Universe(hoppsf,hopdcd)
            group = u.selectAtoms('type *')   
            self.group = group      # group that refers to hopping trajectory
            self.hoptraj = u.dcd    # DCD(!) trajectory object
            self.ts = self.hoptraj.ts
            self.numframes = self.hoptraj.numframes
            self.totaltime = totaltime(self.hoptraj,'ps')
        else:
            raise ValueError('Not sufficient data to create a hopping trajectory.')

    filename = hop.utilities.filename_function

    def next(self):
        """Provides the next time step of a hopping trajectory.

        ts = next()

        If a hopping trajectory file exists then this is
        used. Otherwise, the coordinate trajectory is mapped on the
        fly (which is computationally more expensive).
        """
        if self.hoptraj:
            nextTS = self._read_next_timestep
        else:
            nextTS = self._map_next_timestep
        return nextTS()

    def _map_next_timestep(self):
        """Read next timestep from coordinate trajectory and set up the
        hopping trajectory time step        
        """
        return self._coord2hop(self.traj.next())        

    def _read_next_timestep(self):
        """Read next time step from hopping trajectory"""
        return self.hoptraj.next()

    def write(self,filename,start=None,step=None,delta=None,load=True):
        """Write hopping trajectory as standard dcd file, together with a minimal psf.

        write('hop')

        Arguments:

        load = True     Immediately loads the trajectory so that further
                        calls to next() will use the computed
                        trajectory and don't use expensive mapping.
        
        Ignore the other options and leave them at the
        defaults. Currently, only the whole trajectory is written. For
        visualization one also needs the dummy psf of the group.

        Results:

        filename.dcd and filename.psf

        Note that it is your responsibility to load the hopping
        trajectory and the appropriate psf together as there is very
        limited information stored in the dcd itself.
        """
        set_verbosity(self.verbosity)  # this is stupid
        
        psfname = self.filename(filename,'psf')
        dcdname = self.filename(filename,'dcd')

        # see MDAnalysis/src/dcd/dcd.c for explanations
        if start is None:
            start = self.traj.start_timestep # starting time step for DCD file
        if step is None:
            step = self.traj.skip_timestep   # NSAVC (# ts between written DCD frames)
        if delta is None:
            delta = self.traj.delta          # length of ts (AKMA units)
            
        dcdwriter = MDAnalysis.DCD.DCDWriter(dcdname,self.ts.numatoms,
                                             start,step,delta,
                                             remarks='Hopping trajectory: x=site y=orbit_site z=0')
        for ts in self.map_dcd():
            if ts.frame % 10 == 0:
                msg(3,"Mapping frame %5d/%d  [%5.1f%%]\r" % \
                        (ts.frame,self.numframes,100.0*ts.frame/self.numframes))
            dcdwriter.write_next_timestep(ts)
        dcdwriter.close_trajectory()
        msg(3,"Mapping frame %5d/%d  [%5.1f%%]\r" % \
                (ts.frame,self.numframes,100.0*ts.frame/self.numframes))
        msg(3,'\nFinished writing %s.\n' % dcdname)

        self.write_psf(psfname)
        msg(3,'\nWrote %s.\n' % psfname)

        if load is True:
            self.__init__(filename=filename,verbosity=self.verbosity)
        
    def write_psf(self,filename):
        """Write a dummy psf just for the atoms in the selected group
        so that one can visualize the hopping trajectory.

        write_psf(filename)

        The psf is NOT a fully functional psf. It only contains the
        header and the ATOMS section. It is sufficient to display the
        hopping trajectory in VMD and can be read in by the MDAnalysis
        tools in order to store the atom numbers for the hopping
        trajectory.
        
        ------
        notes
        ------
        Format from psffres.src

        CHEQ:
        II,LSEGID,LRESID,LRES,TYPE(I),IAC(I),CG(I),AMASS(I),IMOVE(I),ECH(I),EHA(I)

        standard format:
          (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8,2G14.6)
          (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8,2G14.6)  XPLOR
        expanded format EXT:
          (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8,2G14.6)
          (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A4,1X,2G14.6,I8,2G14.6) XPLOR
        
        no CHEQ:
        II,LSEGID,LRESID,LRES,TYPE(I),IAC(I),CG(I),AMASS(I),IMOVE(I)

        standard format:
          (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8)
          (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8)  XPLOR
        expanded format EXT:
          (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8)
          (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A4,1X,2G14.6,I8) XPLOR
        """
        # Standard no CHEQ format:
        psf_ATOM_format = '%(iatom)8d %(segid)4s %(resid)-4d %(resname)4s '+\
                          '%(name)-4s %(type)4s %(charge)-14.6f%(mass)-14.4f%(imove)8d\n'
        # This produces:
        #    2114 XWAT 1    TIP3 OH2    75 -0.834000     15.9994              0
        # For some reason, I don't get the same output as I get from Charmm because
        # Charmm centers the value whereas I can only left or right-align them.
        # (This is not a problem as the fields are properly lined up!)
        #    2114 XWAT 1    TIP3 OH2    75  -0.834000       15.9994           0

        psf = open(filename,'w')

        psf.write('PSF\n\n')
        psf.write('%7d !NTITLE\n' % 3)
        psf.write('* Hopping trajectory written by\n'+\
                  '* $Id$\n'+\
                  '* This is NOT a fully functional psf but should work for visualization.\n')
        psf.write('\n')

        psf.write('%6d !NATOM\n' % len(self.tgroup))
        imove = 0    # no fixed atoms
        for atom in self.tgroup:
            # add +1 to atom.number (zero-index but Charmm is 1-indexed) (see PSFParser.py)
            psf.write(psf_ATOM_format % 
                      {'iatom':atom.number+1, 'segid':atom.segid, 'resid':atom.resid,
                       'resname':atom.resname, 'name':atom.name, 'type':atom.type,
                       'charge':atom.charge, 'mass':atom.mass,'imove':imove} )
        # ignore all the other sections (enough for MDAnalysis, VMD, and me)
        psf.close()

    def map_dcd(self,start=None,stop=None,skip=1):
        """Generator to read the trajectory from start to stop and map
        positions to sites.

        ts = map_dcd(**kwargs)

        Arguments:
        start        starting frame number (None means first)
        stop         last frame to read (exclusive) (None means last)
                     (Those are arguments to dcd[start:stop].)
        Iterator Returns:
        ts           hopping trajectory timestep object (iterator)
        """
        # note: iterator + loop is slower than direct loop so I may
        # implement other functions directly with loops and leave the
        # iterator for the user
        if start is not None or stop is not None:
            raise NotImplemented('start/stop do not work on big trajectories')
        if start is None:
            start = 0
        if stop is None:
            stop = self.numframes
            
        self._init_coord2hop()
        #for traj_ts in self.traj[start:stop]:
        for traj_ts in self.traj:              # no slicing for big trajectories
            yield self._coord2hop(traj_ts)    

    def _init_coord2hop(self):
        """Allocate helper arrays for _coord2hop()"""
        # initialization with 'interstitial' is CRUCIAL: throws away first frame
        # and makes sure that we don't keep spurious sites from 1st frame around
        self._sites_last = SITELABEL['interstitial'] * numpy.ones(self.tgroup.numberOfAtoms())
        self._offsites = numpy.empty(self.tgroup.numberOfAtoms(),dtype=bool)
                
    def _coord2hop(self,ts):
        """Translate a single trajectory coordinate frame into a hopping
        trajectory frame and updates the hopping trajectory frame.

        ts          MDAnalysis.dcd.ts  time step object (input coordinate data)

        :Returns:

        hopping ts  Timestep object for the _selected_ atoms with (x=sites y=orbit site z=0)
                    (also updates self.ts so that the HoppingTrajectory instance is uptodate.)
        """
        self.ts.frame = ts.frame   # update the hopping time step        
        coords = numpy.asarray(self.tgroup.coordinates())
        N,D = coords.shape
        
        # Basic nD histograming code from numpy.histogramdd:
        #
        # digitize returns i such that bins[i-1] <= x < bins[i]
        # outliers: i=0 or i=len(bins).
        #    
        # indices[] are NOT map[] indices: to remove the two outlier
        # bins (in the logic of digitize()) we would have to subtract
        # 1 later and also remove indices belonging to outliers. We
        # cheat and add outlier bins to the map (buffered_map[]) and
        # simply label outliers in the trajectory.
        indices = [numpy.digitize(coords[:,i], self.edges[i]) for i in xrange(D)]

        # Using digitize, values that fall on an edge are put in the right bin.
        # For the rightmost bin, we want values equal to the right
        # edge to be counted in the last bin, and not as an outlier.
        for i in xrange(D):
            # Rounding precision
            decimal = int(-numpy.log10(self.dedges[i].min())) +6
            # Find which points are on the rightmost edge.
            on_edge = numpy.where(
                numpy.around(coords[:,i], decimal) == \
                numpy.around(self.edges[i][-1], decimal))[0]
            # Shift these points one bin to the left.
            indices[i][on_edge] -= 1
            
        # indices contains the outliers at index 0 and len(edges[i])
        # To make things simpler, we expand the map with a outlier zone,
        # label outliers with -1 and then index into the buffered_map

        # fill the core of the buffered map
        core = D*[slice(1,-1)]
        self.buffered_map[core] = self.map          # not very expensive

        # Note that now indices[m] corresponds to group[m]
        # (It's important that the order is preserved)
        # The group is not needed right now, though.
        #
        # pos[:,0] = site(t), pos[:,1] = orbit site, pos[:,2] = 0 (unused)
        pos = self.ts._pos     # assign slices to avoid loop (thanks to Naveen)
        pos[:,0] = [self.buffered_map[indices[0][iatom],indices[1][iatom],indices[2][iatom]] \
                    for iatom in xrange(N)]
        s = pos[:,0]
        self._offsites[:] = (s == SITELABEL['interstitial']) | (s == SITELABEL['outlier'])
        pos[:,1] = s      # particles in interstital and outliers are assigned their previous site
        pos[self._offsites,1] = self._sites_last[self._offsites]
        pos[:,2] = 0

        # _sites_last[] was initialized to 'interstitial': this ensures proper accounting 
        # for all later steps (because 'interstitial' is thrown away at the analysis stage)
        self._sites_last[:] = pos[:,1]  # save orbit sites for next step
        return self.ts

    def iterator(self):
        return self.__iter__()

    def __iter__(self):
        if self.hoptraj:
            for ts in self.hoptraj:
                yield ts
        else:
            self._init_coord2hop()
            for traj_ts in self.traj:
                yield self._coord2hop(traj_ts)    

class TAPtrajectory(object):
    """Provides a Time-Averaged Position (TAP) version of the input
    trajectory. The method is described in Henchman and McCammon, J
    Comp Chem 23 (2002), 861 doi:10.1002/jcc.10074

    :Attributes:

    ts               MDAnalysis.Timestep object
    numframes        number of frames in TAP trajectory
    group            AtomGroup of atoms that are tracked


    :Methods:

    ## [start:stop]         object can be used as an iterator over the
    ##                      hopping trajectory (disabled due to dcdreader bug)
    next()      advances time step in the hopping trajectory
    map_dcd()            iterator that updates the ts and maps the trajectory
                         coordinates to site labels
    _map_next_timestep()  map next coordinate trajectory step to hopping time step
    _read_next_timestep() read next timestep from hopping trajectory

    write()              write the hopping trajectory to a dcd file + psf
    """

    def __init__(self,trajectory=None,group=None,TAPradius=2.8,TAPsteps=3,
                 filename=None,dcd=None,psf=None,fixtrajectory=None,verbosity=3):
        """A TAP trajectory object converts a trajectory into a TAP trajectory.

        Create from a coordinate trajectory of a group of water residues:

          u = MDAnalysis.Universe(psf,dcd)
          water = u.selectAtoms('resname TIP*')  # see NOTE below!!
          water = u.selectAtoms('name OH2')      # better, see NOTE below!!
          h = TAPtrajectory(trajectory=u.dcd,group=water)

        Load from a saved hopping trajectory (in dcd format with dummy psf)

          h = TAPtrajectory(dcd='TAP.dcd',psf='TAP.psf')

        The given atom group is filtered according to the Time-Averaged Positon
        algorithm (Henchman and McCammon, J Comp Chem 23 (2002), 861). Original
        positions are replaced by their TAPs: A particles last position (TAP)
        is retained unless it has moved farther than TAPradius from its TAP
        measured by its root mean square distance over the last TAPsteps
        frames.

        One can use a TAP filtered trajectory 'on-the-fly' to build the density:
          
          u = Universe(psf,dcd)
          oxy = u.selectAtoms('name OH2')
          TAP = TAPtrajectory(u.dcd,oxy)
          u.dcd = TAP.dcd    # <--- replace orig dcd with TAP !!
          dens = hop.sitemap.density_from_Universe(u,atomselection='name OH2')

        NOTE: In the current implementation residues are often ripped apart
        because all coordinates are processed independently. It is recommended
        to only do TAP on the water oxygens (for speed). This will create a
        trajectory in which hydrogens are always ripped from the oxygen but
        this trajectory is ONLY being used for creating a density from those
        oxygen using hop.sitemap.build_density().

        (This could be fixed at the cost of speed; in this case TAP would be done
        on the centre of mass and the whole residue would be translated.)

        :Arguments:

        trajectory       MDAnalysis.dcd trajectory instance
        group            MDAnalysis.group instance (from the same Universe as trajectory)
        TAPradius        particles are considered to be on the TAP as long as they 
                         haven't moved farther than TAPradius over the last TAPsteps frames
        TAPsteps         RMS distance of particle from TAP over TAPsteps is compared
                         to TAPradius
        dcd              dcd written by write()
        psf              psf written by write() (or write_psf())
        filename         or simply provide one filename prefix for psf and dcd

        fixtrajectory    dictionary with attributes of a dcd object and new
                         values; used to provide correct values after using
                         a catdcd-generated trajectory (hack!), e.g.
                         fixtrajectory = {'delta':10.22741474887299}
        
        verbosity        show status messages for >= 3
        """
        self.verbosity = verbosity
        set_verbosity(self.verbosity)

        if not (trajectory is None or group is None):
            self.traj  = trajectory                      # MDAnalysis.Universe.trajectory
            self.tgroup = group                          # atom selection for trajectory
            self.tgroup_indices = self.tgroup.indices()  # cache indices
            if not isinstance(self.tgroup,MDAnalysis.AtomGroup.AtomGroup):
                raise TypeError('group must be a <AtomGroup>, eg MDAnalyis.Universe.selectAtoms().')
            self.universe = self.tgroup.atoms[0].universe  # Universe of dcd and group (hackish..)
            if isinstance(fixtrajectory,dict):
                for attr,val in fixtrajectory.items():
                    if not hasattr(trajectory,attr):
                        raise AttributeError('fixtrajectory: dcd object does not have attribute "'\
                                             +str(attr)+'"')
                    trajectory.__dict__[attr] = val 
            self.totaltime = totaltime(trajectory,'ps')
            self.traj.rewind()             # make sure to start from frame 0
            self.ts = self.traj.ts         # output will look like input (no copy, see _coord2TAP!)
            self.TAPtraj = None            # no TAP trajectory available
            self.TAPradius = TAPradius
            self.TAPsteps = TAPsteps
            # store last TAPsteps in __lastframes
            self.__lastframes = hop.utilities.Ringbuffer(self.TAPsteps)
            # store the last TAP coordinates: initialized here
            self.__currentTAP = self.tgroup.coordinates().copy()
            # fake DCD object that can be slotted into another universe
            self.dcd_attributes = {}
            for k in ['dcdfilename','delta','filename','fixed','numframes',
                      'numatoms','periodic','remarks',
                      'skip','skip_timestep','start_timestep']:
                self.dcd_attributes[k] = self.traj.__dict__[k]
            self.dcd = ThinDCDReader(self)
            self.numframes = self.dcd_attributes['numframes']
        elif not (dcd is None or psf is None) or filename is not None:
            # read from dcd
            try:
                self.traj,self.tgroup
            except AttributeError:
                self.traj,self.tgroup = [None] * 2
            if filename is not None:
                psf = self.filename(filename,'psf')
                dcd = self.filename(filename,'dcd')
            u = MDAnalysis.Universe(psf,dcd)
            group = u.selectAtoms('type *')   # TODO: why do I need this?
            self.group = group      # group that refers to hopping trajectory
            self.TAPtraj = u.trajectory    # DCD trajectory object
            self.ts = self.TAPtraj.ts
            self.numframes = self.TAPtraj.numframes
            self.totaltime = totaltime(self.TAPtraj,'ps')
            # DCD object that can be slotted into another universe
            self.dcd = u.trajectory
        else:
            raise ValueError('Not sufficient data to create a TAP trajectory.')

    filename = hop.utilities.filename_function

    def next(self):
        """Provides the next time step of a TAP trajectory.

        ts = next()

        If a TAP trajectory file exists then this is used. Otherwise,
        the coordinate trajectory is mapped on the fly (which is
        computationally more expensive).
        """
        if self.TAPtraj:
            nextTS = self._read_next_timestep
        else:
            nextTS = self._map_next_timestep
        return nextTS()

    def rewind(self):
        if self.TAPtraj:
            self.TAPtraj.rewind()
        else:
            self.traj.rewind()

    def _map_next_timestep(self):
        """Read next timestep from coordinate trajectory and set up the
        TAP trajectory time step        
        """
        return self._coord2TAP(self.traj.next())        

    def _read_next_timestep(self):
        """Read next time step from a TAP trajectory on disk"""
        return self.TAPtraj.next()

    def write(self,filename,start=None,step=None,delta=None,load=True):
        """Write hopping trajectory as standard dcd file.

        write('TAP')

        :Arguments:

        load = True     Immediately loads the trajectory so that further
                        calls to next() will use the computed
                        trajectory and don't use expensive mapping.
        
        Ignore the other options and leave them at the defaults. Currently,
        only the whole trajectory is written. All atoms in the original
        trajectory are written to the output so you should be able to use your
        original psf file.

        NOTE: Fixed atoms are possibly not accounted for properly.

        Note that it is your responsibility to load the TAP trajectory and the
        appropriate psf together as there is very limited information stored in
        the dcd itself.
        """
        set_verbosity(self.verbosity)  # this is stupid
        
        psfname = self.filename(filename,'psf')
        dcdname = self.filename(filename,'dcd')

        # see MDAnalysis/src/dcd/dcd.c for explanations
        if start is None:
            start = self.traj.start_timestep # starting time step for DCD file
        if step is None:
            step = self.traj.skip_timestep   # NSAVC (# ts between written DCD frames)
        if delta is None:
            delta = self.traj.delta          # length of ts (AKMA units)
            
        dcdwriter = MDAnalysis.DCD.DCDWriter(dcdname,self.ts.numatoms,
                                             start,step,delta,
                                             remarks='TAP trajectory')
        for ts in self.map_dcd():
            if ts.frame % 10 == 0:
                msg(3,"Mapping frame %5d/%d  [%5.1f%%]\r" % \
                        (ts.frame,self.numframes,100.0*ts.frame/self.numframes))
            dcdwriter.write_next_timestep(ts)
        dcdwriter.close_trajectory()
        msg(3,"Mapping frame %5d/%d  [%5.1f%%]\r" % \
                (ts.frame,self.numframes,100.0*ts.frame/self.numframes))
        msg(3,'\nFinished writing %s.\n' % dcdname)

        if load is True:
            self.TAPtraj = MDAnalysis.DCD.DCDReader(dcdname)
            self.dcd = self.TAPtraj
        
    def map_dcd(self,start=None,stop=None,skip=1):
        """Generator to read the trajectory from start to stop and map
        positions to TAP sites.

        ts = map_dcd(**kwargs)

        Arguments:
        start        starting frame number (None means first)
        stop         last frame to read (exclusive) (None means last)
                     (Those are arguments to dcd[start:stop].)
        Iterator Returns:
        ts           hopping trajectory timestep object (iterator)
        """
        # note: iterator + loop is slower than direct loop so I may
        # implement other functions directly with loops and leave the
        # iterator for the user
        if start is not None or stop is not None:
            raise NotImplemented('start/stop do not work on big trajectories')
        if start is None:
            start = 0
        if stop is None:
            stop = self.numframes

        #for traj_ts in self.traj[start:stop]:
        for traj_ts in self.traj:              # no slicing for big trajectories(ERRORS!)
            yield self._coord2TAP(traj_ts)    
                
    def _coord2TAP(self,ts):
        """Translate a single trajectory coordinate frame into a TAP
        trajectory frame and update the TAP trajectory frame.

        ts        MDAnalysis.dcd.ts  time step object

        Only the selection's coordinates are TAP-filtered.
        :Returns:

        hopping  ts
        """
        # Modify the original frame in place and use it as the new frame. This
        # should work in most instances unless one wants to immediately compare old
        # and new frame.
        self.ts = ts          # update the TAP time step 
        
        # only work on the selected coordinates (memory efficiency but
        # slower??)  (I didn't manage to always work on a reference to the
        # coords; this would avoid having to patch back the altered coordinates
        # into the whole coord set, see below.)
        coords = self.tgroup.coordinates()      # makes a new copy
        self.__lastframes.append(coords.copy()) # remember last TAPsteps frames

        # calculated RMS distance for last TAPsteps from current TAP for all coords
        # (current coords are part of the running average over __lastframes)
        # __currentTAP is initialized in __init__ to the first frame
        d = numpy.average(
            numpy.sqrt(
                numpy.sum((numpy.asarray(self.__lastframes) - self.__currentTAP), axis=0)**2),
            axis=1)
        onTAP = (d <= self.TAPradius)              # particles that did not move far    
        coords[onTAP] = self.__currentTAP[onTAP]   # reset to TAP
        self.__currentTAP[:] = coords              # remember TAP for next frame
        # patch  selected coordinates back into full coordinate set
        #    u.dcd.ts._pos[w.indices()] = new_coord_array    # WORKS (no copy)
        #    x = u.dcd.ts._pos[w.indices()]                  # FAILS (copy involved)
        #    x[:] = new_coord_array[:]                       # 
        self.ts._pos[self.tgroup_indices] = self.__currentTAP
        return self.ts

    def iterator(self):
        return self.__iter__()

    def __iter__(self):
        if self.TAPtraj:
            for ts in self.TAPtraj:
                yield ts
        else:
            for traj_ts in self.traj:
                yield self._coord2TAP(traj_ts)    




def RMS_fit_trj(traj,ref,select='backbone',filename=None,prefix='rmsfit_',verbosity=3):
    """RMS-fit trajectory to a reference structure using a selection.

    RMS_fit_trj(traj, ref, 'backbone or name CB or name OT*')

    :Input:
    traj       trajectory, MDAnalysis universe object
    ref        reference coordinates; MDAnalysis universe object
               (uses the current time step of the object)
    filename   file name for the RMS-fitted trajectory or pdb; defaults to the 
               original dcd filename (from traj) with prefix prepended
    prefix     prefix for autogenerating the new dcd filename
    select     any valid selection string for MDAnalysis.AtomGroup.selectAtom
               that produces identical selections in traj and ref
               or dictionary {'reference':sel1, 'target':sel2}.
               The fasta2select() function returns such a dictionary based on a
               clustalW or STAMP sequence alignment.
    verbosity  at level 3, display progress; at 5, also computes new/old RMSD

    Both reference and trajectory must be MDAnalysis.Universe
    instances. If they contain a dcd then this is used. Alternatively,
    they use its pdb object; if neither dcd nor pdb exist an
    AttributeError is raised.  If traj contains a dcd then a dcd is
    written, otherwise only a single pdb is written.

    Note that this function is flexible but slow. In order to fit long
    trajectory it is suggested to only align one frame and then use
    that frame with a more efficient program (that otherwise could not
    easily deal with aligning different molecules).
    """

    # Charmm (using OB's rmsfit.inp script):
    # charmm < /Users/oliver/Biop/Library/Charmm/analysis/trajectory/rmsfit.inp \
    #   HEADER=str/header.str PSF=inp/crbp_apo.psf DCD=trj/1opa_100ps_ppc.dcd \
    #   DCD_RMS=trj/rmsfit2_1opa_100ps_ppc.dcd  RECENTER=0 \
    #   REF_PDB=coord/rmsfit_1opa_a_xtal.pdb 

    import numpy
    import os.path

    set_verbosity(verbosity)

    # setup new dcd
    try:
        frames = traj.dcd
    except AttributeError:
        try:
            frames = traj.pdb
        except AttributeError:
            raise AttributeError('traj must either contain a DCD or a PDB.')

    if filename is None:
        path,fn = os.path.split(frames.filename)
        filename = os.path.join(path,prefix+fn)
    if type(select) is str:
        select = {'reference':select,'target':select}

    if isinstance(frames,MDAnalysis.DCD.DCDReader):
        writer = MDAnalysis.DCD.DCDWriter(filename,frames.numatoms,
                                          frames.start_timestep,
                                          frames.skip_timestep,
                                          frames.delta,
                                          remarks='RMS fitted trajectory to ref')
    elif isinstance(frames,MDAnalysis.PDB.PDBReader):
        writer = MDAnalysis.PDB.PDBWriter(filename, universe=traj,
                                          remarks='RMS fitted pdb frame to ref')
    else:
        raise TypeError('traj must either contain a DCD or a PDB.')

    ref_atoms = ref.selectAtoms(select['reference'])
    traj_atoms = traj.selectAtoms(select['target'])
    if len(ref_atoms) != len(traj_atoms):
        raise SelectionError("Reference and trajectory atom selections do not contain "+
                             "the same number of atoms: N_ref=%d, N_traj=%d" % \
                             (len(ref_atoms), len(traj_atoms)))
    msg(4,"RMS-fitting on %d atoms.\n" % len(ref_atoms))
    mass_mismatches = (numpy.absolute(ref_atoms.masses() - traj_atoms.masses()) > 0)
    if numpy.any(mass_mismatches):
        # diagnostic output:
        msg(0,"Atoms: reference | trajectory\n")
        for ar,at in zip(ref_atoms,traj_atoms):
            if ar.name != at.name:
                msg(0,"%4s %3d %3s %3s %6.3f  |  %4s %3d %3s %3s %6.3f\n" %  \
                      (ar.segid, ar.resid, ar.resname, ar.name, ar.mass,
                       at.segid, at.resid, at.resname, at.name, at.mass,))
        raise SelectionError("Inconsistent selections, masses don't match; "
                             "mis-matching atoms are shown above.")
    del mass_mismatches
    masses = ref_atoms.masses()
    
    # reference centre of mass system
    # (compatibility with pre 1.0 numpy: explicitly cast coords to float32)
    ref_com = ref_atoms.centerOfMass().astype(numpy.float32)
    ref_coordinates = ref_atoms.coordinates() - ref_com

    # allocate the array for selection atom coords
    traj_coordinates = traj_atoms.coordinates().copy()

    # R: rotation matrix that aligns r-r_com, x~-x~com   
    #    (x~: selected coordinates, x: all coordinates)
    # Final transformed traj coordinates: x' = (x-x~_com)*R + ref_com
    for ts in frames:
        # shift coordinates for rotation fitting
        # selection is updated with the time frame
        x_com = traj_atoms.centerOfMass().astype(numpy.float32)
        traj_coordinates[:] = traj_atoms.coordinates() - x_com
        R = numpy.matrix(MDAnalysis.rms_fitting.rms_rotation_matrix(traj_coordinates,ref_coordinates,masses),dtype=numpy.float32)
        # Transform each atom in the trajectory (use inplace ops to avoid copying arrays)
        ts._pos   -= x_com
        ts._pos[:] = ts._pos * R # R acts to the left & is broadcasted N times.
        ts._pos   += ref_com
        writer.write_next_timestep(ts)
        if msg(5):
            rmsd_old = rmsd(ref_atoms.coordinates(),traj_coordinates)
            rmsd_new = rmsd(ref_atoms.coordinates(),traj_atoms.coordinates())
            msg(5,"Fitted frame %5d/%d  [%5.1f%%]  %5.2fA --> %5.2fA  |translation|=%.2fA\r" % \
                    (ts.frame,frames.numframes,100.0*ts.frame/frames.numframes,
                     rmsd_old, rmsd_new, rmsd(x_com,ref_com)) )
        else:
            if ts.frame % 10 == 0:
                msg(3,"Fitted frame %5d/%d  [%5.1f%%]\r" % \
                        (ts.frame,frames.numframes,100.0*ts.frame/frames.numframes))
    msg(3,"Fitted frame %5d/%d  [%5.1f%%]\r" % \
            (ts.frame,frames.numframes,100.0*ts.frame/frames.numframes))
    msg(3,"\nWrote RMS-fitted coordinates to file '%s'\n" % filename)

def rmsd(a,b):
    """Returns RMSD between two coordinate sets a and b."""
    return numpy.sqrt(numpy.sum(numpy.power(a-b,2))/a.shape[0])

def fasta2select(fastafilename,is_aligned=False,
                 ref_resids=None, target_resids=None,
                 ref_offset=0,target_offset=0,verbosity=3):
    """Align two sequences from a fasta file and construct a MDAnalysis
    selection string of the common atoms.

    select_dict = fasta2select(fastafilename)

    fastafilename contains the the two un-aligned sequences in FASTA
    format. The reference is assumed to be be the first sequence, the
    target the second. Clustalw produces a pairwise alignment (which
    is written to a file with suffix .aln).  The output contains atom
    selection strings that select the same atoms in the two
    structures.

    Unless ref_offset and/or target_offset are specified, the resids
    in the structure are assumed to correspond to the positions in the
    un-aligned sequence, namely the first residue has resid == 1.

    In more complicated cases (e.g. when the resid numbering in the
    structure/psf has gaps due to missing parts), simply provide the
    sequenxe of resids as they appear in the psf in ref_resids or
    target_resids, eg

       target_resids = [a.resid for a in trj.selectAtoms('name CA')]
 
    (This translation table IS combined with any value for offset!)

    :Input:
    fastafilename    FASTA file with first sequence as reference and
                     second the one to be aligned (ORDER IS IMPORTANT!)
    is_aligned       False: run clustalw for sequence alignment; True: use
                     the alignment in the file (e.g. from STAMP)
    ref_offset       add this number to the column number in the FASTA
                     to get the original residue number
    target_offset    same for the target
    ref_resids       sequence of resids as they appear in the
                     reference structure
    target_resids    sequence of resids as they appear in the target

    :Output:
    select_dict      dictionary with 'reference' and 'target' selection
                     string that can be used immediately in 
                     RMS_fit(...,select=select_dict).
    """
    import Bio
    import numpy

    set_verbosity(verbosity)
    
    if is_aligned:
        import Bio.SeqIO, Bio.Alphabet
        protein_gapped = Bio.Alphabet.Gapped(Bio.Alphabet.IUPAC.protein)
        fasta = open(fastafilename)
        try:
            alignment = Bio.SeqIO.to_alignment(
                Bio.SeqIO.FastaIO.FastaIterator(fasta,alphabet=protein_gapped),
                alphabet=protein_gapped)
        finally:
            fasta.close()
        msg(3,"Using provided alignment, %s\n" % fastafilename)
    else:
        import Bio.Clustalw
        import os.path
        filepath,ext = os.path.splitext(fastafilename)
        alnfilename = filepath + '.aln'
        cline = Bio.Clustalw.MultipleAlignCL(fastafilename)
        cline.set_output(alnfilename)
        cline.set_type('protein')
        alignment = Bio.Clustalw.do_alignment(cline)
        msg(3,"Using clustalw sequence alignment, %s.\n" % alnfilename)

    nseq = len(alignment._records)    # the stupid class should provide __len__ !
    if nseq != 2:
        raise ValueError("Only two sequences in the alignment can be processed.")

    orig_resids = [ref_resids,target_resids] # implict assertion that
                                             # we only have to sequenceses in the alignment
    offsets = [ref_offset,target_offset]
    for iseq,a in enumerate(alignment):      # need iseq index to change orig_resids
        if orig_resids[iseq] is None:
            # build default: assume consecutive numbering of all
            # residues in the alignment
            GAP = a.seq.alphabet.gap_char
            length = len(a.seq) - a.seq.count(GAP)
            orig_resids[iseq] = numpy.arange(1,length+1)
        else:
            orig_resids[iseq] = numpy.asarray(orig_resids[iseq])
    # add offsets to the sequence <--> resid translation table
    seq2resids = [resids + offset for resids,offset in zip(orig_resids,offsets)]
    del orig_resids
    del offsets

    def resid_factory(alignment,seq2resids):
        """Return a function that gives the resid for a position ipos in
        the nseq'th alignment.
        
        resid = resid_factory(alignment,seq2resids)
        r = resid(nseq,ipos)

        It is based on a look up table that translates position in the
        alignment to the residue number in the original
        sequence/structure.
        
        The first index of resid() is the alignmment number, the
        second the position in the alignment.

        seq2resids translates the residues in the sequence to resid
        numbers in the psf. In the simplest case this is a linear map
        but if whole parts such as loops are ommitted from the protein
        the seq2resids may have big gaps.  

        Format: a tuple of two numpy arrays; the first array is for
        the reference, the second for the target, The index in each
        array gives the consecutive number of the amino acid in the
        sequence, the value the resid in the structure/psf.

        Note: assumes that alignments have same length and are padded if
        necessary.
        """
        # could maybe use Bio.PDB.StructureAlignment instead?
        nseq = len(alignment._records)
        t = numpy.zeros((nseq,alignment.get_alignment_length()),dtype=int)
        for iseq,a in enumerate(alignment):
            GAP = a.seq.alphabet.gap_char
            t[iseq,:] = seq2resids[iseq][numpy.cumsum(numpy.where(
                        numpy.array(list(a.seq))==GAP,0,1)) - 1]  
            # -1 because seq2resid is index-1 based (resids start at 1)

        def resid(nseq,ipos):
            return t[nseq,ipos]
        return resid
    resid = resid_factory(alignment,seq2resids)

    res_list = []     # collect individual selection string
    # could collect just resid and type (with/without CB) and
    # then post-process and use ranges for continuous stretches, eg
    # ( resid 1:35 and ( backbone or name CB ) ) or ( resid 36 and backbone ) ...

    GAP = alignment.get_seq_by_num(0).alphabet.gap_char # should be same for all seqs
    for ipos in xrange(alignment.get_alignment_length()):
        aligned = list(alignment.get_column(ipos))
        if GAP in aligned:
            continue       # skip residue
        template = "resid %i"
        if 'G' not in aligned:
            # can use CB
            template += " and ( backbone or name CB )"
        else:
            template += " and backbone"
        template = "( "+template+" )"

        res_list.append([template % resid(iseq,ipos) for iseq in xrange(nseq)])

    sel = numpy.array(res_list).transpose()

    ref_selection =  " or ".join(sel[0])
    target_selection =  " or ".join(sel[1])
    return {'reference':ref_selection, 'target':target_selection}

class ThinDCDReader(MDAnalysis.DCD.DCDReader):
    """DCD-like object that supports a subsection of the DCDReader
    interface such as iteration over frames and most attributes. The
    important part is that the __iter__() method is overriden to
    provide data from another source. This allows a filter architecture
    for trajectories."""
    # Right now specifically designed for TAPtrajectory class.
    def __init__(self,datafeeder):
        D = datafeeder
        # should have as attributes:
        # ['dcdfilename','delta','filename','fixed','numframes', 'numatoms','periodic',
        #  'remarks', 'skip','skip_timestep','start_timestep']
        self.__dict__.update(D.dcd_attributes)
        self.dcdfile = None  # no file is linked; with None, __del__ will be happy
        # use the classes/methods from the feeder class:
        self.ts = D.ts
        self.__iter__ = D.__iter__        
        self.next = D.next  # feeder needs next()
        self.rewind = D.rewind
    def __getitem__(self,frame):
        """Slow sequential 'read forward' implementation."""
        for ts in self:
            if ts.frame == frame+1:  # frames are 1-based
                break
        return ts
    def rewind(self):
        self[0]
    def timeseries(self,*args,**kwargs):
        raise NotImplementedError
    def correl(self,*args,**kwargs):
        raise NotImplementedError
    def close_trajectory(self):
        pass
    
