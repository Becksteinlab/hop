# Hop --- a framework to analyze solvation dynamics from MD simulations
# Copyright (c) 2007-2010 Oliver Beckstein <orbeckst@gmail.com>
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

"""
Generating densities from trajectories --- :mod:`hop.density`
=============================================================

As an input a trajectory is required that

1. Has been centered on the protein of interest.
2. Has all molecules made whole that have been broken across periodic
   boundaries.
3. Has the solvent molecules remap so that they are closest to the
   solute (this is important when using funky unit cells such as
   dodechedra or truncated octahedra).

Classes and functions
---------------------
"""
from __future__ import absolute_import

import sys
import os, os.path
import errno
import cPickle
import warnings

import numpy
import MDAnalysis
from gridData import OpenDX

from . import MissingDataError, InconsistentDataWarning
from . import constants
from .utilities import msg,set_verbosity,get_verbosity, flatten, sorted, \
     DefaultDict, fixedwidth_bins, iterable, asiterable, CustomProgressMeter
from .sitemap import Density


import logging
logger = logging.getLogger("hop.density")


class DensityCollector(object):
    """Collect subsequent coordinate frames to build up a :class:`Density`."""
    use_kdtree = True
    def __init__(self, name, universe, **kwargs):
        self.name = name
        try:
            universe.selectAtoms('all')
            universe.trajectory.ts
        except AttributeError:
            errmsg = "DensityCollector: The universe must be a proper MDAnalysis.Universe instance."
            logger.fatal(errmsg)
            raise TypeError(errmsg)
        self.universe = u = universe
        self.delta = kwargs.pop('delta', 1.0)
        self.atomselection = kwargs.pop('atomselection', 'name OH2')
        self.cutoff = kwargs.pop('cutoff', 3.5)
        self.soluteselection = kwargs.pop('soluteselection', None) #'protein and not name H*')
        self.padding = kwargs.pop('padding', 2.0)
        self.metadata = kwargs.pop('metadata', {})
        self.parameters = kwargs.pop('parameters',{})  # for advanced fiddling...
        # define the self.current_coordinates() function ... monkey patching!
        if self.cutoff > 0 and self.soluteselection is not None:
            # special fast selection for '<atomsel> not within <cutoff> of <solutesel>'
            notwithin_coordinates = notwithin_coordinates_factory(u,self.atomselection,self.soluteselection,
                                                                  self.cutoff,use_kdtree=self.use_kdtree)
            self.current_coordinates = notwithin_coordinates
            self.mode = "BULK"
        else:
            group = u.selectAtoms(self.atomselection)
            self.current_coordinates = group.coordinates
            self.mode = "SOLVENT"
        coord = self.current_coordinates()
        logger.info("%-10s: Selected %d atoms out of %d atoms (%s) from %d total.",
                    self.name, coord.shape[0],len(u.selectAtoms(self.atomselection)),
                    self.atomselection,len(u.atoms))

        self.__complete = False

    def init_histogram(self, **kwargs):
        # needs to be done separately because we might need additional information
        # after init (at least I cannot think of a better way...)
        smin = kwargs.pop("smin", self.min_coordinates(padding=self.padding))
        smax = kwargs.pop("smax", self.max_coordinates(padding=self.padding))

        BINS = fixedwidth_bins(self.delta, smin, smax)
        self.arange = zip(BINS['min'],BINS['max'])
        self.bins = BINS['Nbins']

        # create empty grid with the right dimensions (and get the edges)
        grid,edges = numpy.histogramdd(numpy.zeros((1,3)), bins=self.bins,
                                       range=self.arange, normed=False)
        grid *= 0.0
        h = grid.copy()

        self.grid = grid
        self.edges = edges
        self._h = h         # temporary for accumulation

    def min_coordinates(self, **kwargs):
        return numpy.min(self.current_coordinates(), axis=0) - kwargs.pop('padding', self.padding)

    def max_coordinates(self, **kwargs):
        return numpy.max(self.current_coordinates(), axis=0) + kwargs.pop('padding', self.padding)

    def collect(self):
        assert hasattr(self, 'grid'), "init_histogram() must be called first"
        coord = self.current_coordinates()
        if len(coord) > 0:
            self._h[:],self.edges[:] = numpy.histogramdd(coord, bins=self.bins, range=self.arange, normed=False)
            self.grid += self._h  # accumulate average histogram
        return len(coord)

    def finish(self):
        if self.isComplete():
            return
        u = self.universe
        numframes = u.trajectory.numframes / u.trajectory.skip
        self.grid /= float(numframes)
        self.__complete = True

    def Density(self):
        if not hasattr(self, 'grid'):
            errmsg = "DensityCollector.Density(): No data for density available. Run collect() first."
            logger.error(errmsg)
            raise MissingDataError(errmsg)
        u = self.universe
        metadata = self.metadata
        metadata['collector'] = self.name
        metadata['collector_mode'] = self.mode
        metadata['psf'] = u.filename             # named psf for historical reasons: any topol
        metadata['dcd'] = u.trajectory.filename  # named dcd for historical reasons: any traj
        metadata['atomselection'] = self.atomselection
        metadata['numframes'] = u.trajectory.numframes
        metadata['dt'] = u.trajectory.dt    # in ps for default MDAnalysis
        # totaltime should be in MDAnalysis!
        metadata['totaltime'] = round(u.trajectory.numframes * metadata['dt'] * u.trajectory.skip_timestep, 3)
        metadata['time_unit'] = MDAnalysis.core.flags['time_unit']  # just to make sure we know it...
        metadata['dcd_skip'] = u.trajectory.skip_timestep  # frames
        metadata['dcd_delta'] = u.trajectory.delta         # in native units (?)
        if self.mode == 'BULK':
            metadata['soluteselection'] = self.soluteselection
            metadata['cutoff'] = self.cutoff             # in Angstrom

        parameters = self.parameters
        parameters['isDensity'] = False             # must override

        # Density automatically converts histogram to density for isDensity=False
        g = Density(grid=self.grid, edges=self.edges,
                    unit=dict(length=MDAnalysis.core.flags['length_unit']),
                    parameters=parameters, metadata=metadata)
        logger.info("%-10s: Histogram completed (initial density in %s**-3)",
                    self.name, MDAnalysis.core.flags['length_unit'])
        return g

    def isComplete(self):
        return self.__complete

    def __repr__(self):
        if self.mode == "BULK":
            return "<DensityCollector %(name)r, delta=%(delta).1f A: "\
                "'%(atomselection)s and not around %(cutoff).1f (%(soluteselection)s)'>" % vars(self)
        else:
            return "<DensityCollector %(name)r, delta=%(delta).1f A: %(atomselection)r>" % vars(self)


class DensityCreator(object):
    modes = ("all", "bulk", "solvent")
    defaults = {'cutoff': 3.5,
                'soluteselection': "protein and not name H*",
                'delta':1.0, 'atomselection': "name OH2",
                'padding': 2.0,
                }

    def __init__(self, *args, **kwargs):
        """Create a density grid from a trajectory.

           density_from_trajectory(PSF, DCD, delta=1.0, atomselection='name OH2', ...) --> density

        or

           density_from_trajectory(PDB, XTC, delta=1.0, atomselection='name OH2', ...) --> density

        :Arguments:
          psf/pdb/gro
                topology file
          dcd/xtc/trr/pdb
                trajectory; if reading a single PDB file it is sufficient to just provide it
                once as a single argument

        :Keywords:
          mode
                'solvent', 'bulk' or 'all' ('all' does both 'solvent' and \bulk' at the
                same time and thus :meth:`DensityCreator.Density`` returns a list of
                densities; this saves time!)  ['all']
          atomselection
                selection string (MDAnalysis syntax) for the species to be analyzed
                ["name OH2"]
          delta
                approximate bin size for the density grid in Angstroem (same in x,y,z)
                (It is slightly adjusted when the box length is not an integer multiple
                of delta.) [1.0]
          metadata
                dictionary of additional data to be saved with the object
          padding
                increase histogram dimensions by padding (on top of initial box size)
                in Angstroem [2.0]
          soluteselection
                MDAnalysis selection for the solute, e.g. "protein" [``None``]
          cutoff
                With *cutoff*, select '<atomsel> NOT WITHIN <cutoff> OF <soluteselection>'
                (Special routines that are faster than the standard AROUND selection) [0]
          verbosity: int
                level of chattiness; 0 is silent, 3 is verbose [3]

        :Returns: :class:`hop.sitemap.Density`

        :TODO:
          * Should be able to also set skip and start/stop for data collection.

        .. Note::
            * In order to calculate the bulk density, use

                  atomselection='name OH2',soluteselection='protein and not name H*',cutoff=3.5

              This will select water oxygens not within 3.5 A of the protein heavy atoms.
              Alternatively, use the VMD-based  :func:`density_from_volmap` function.
            * The histogramming grid is determined by the initial frames min and max.
            * metadata will be populated with psf, dcd, and a few other items.
              This allows more compact downstream processing.

        """
        _kwargs = self.defaults.copy()
        _kwargs.update(kwargs)
        kwargs = _kwargs
        # workaround for python 2.5 *args,**kwargs only allowed:
        universe_kwargs = {'permissive':kwargs.pop('permissive',False)}
        self.universe = MDAnalysis.asUniverse(*args, **universe_kwargs)
        self.mode = kwargs.pop("mode", "all")   # 'all' runs modes[1:]
        if not self.mode in self.modes:
            errmsg = "DensityCreator: mode must be one of %r, not %r" % (self.modes, self.mode)
            logger.fatal(errmsg)
            raise ValueError(errmsg)
        if self.mode == "all":
            modes = self.modes[1:]
        else:
            modes = [self.mode]
        self.collectors = []
        min_coords = []
        max_coords = []
        for mode in modes:
            modeargs = kwargs.copy()
            if mode == "solvent":
                modeargs['soluteselection'] = None
                modeargs['cutoff'] = 0
            c = DensityCollector(mode, self.universe, **modeargs)
            self.collectors.append(c)
            min_coords.append(c.min_coordinates())   # with default padding from modeargs
            max_coords.append(c.max_coordinates())
        # determine maximum bounding box from initial positions of solvent
        # (add generous padding... probably more than my default 2 A)
        smin = numpy.sort(min_coords, axis=0)[0]   # the three smallest values
        smax = numpy.sort(max_coords, axis=0)[-1]  # the three largest values
        for c in self.collectors:
            c.init_histogram(smin=smin, smax=smax)  # also guarantees compatible grid
        self.densities = {}   # densities will be stored with mode as key

    def create(self):
        """Loop through trajectory and build all densities.

        .. SeeAlso::

           :class:`DensityCollector`
        """
        u = self.universe
        pm = CustomProgressMeter(u.trajectory.numframes, interval=10,
                                 format="Histogramming %(other)s atoms in frame %(step)5d/%(numsteps)d  [%(percentage)5.1f%%]\r")

        for ts in u.trajectory:
            status = ""
            for c in self.collectors:
                natoms = c.collect()
                status += ("%s=%d " % (c.name, natoms))
            pm.echo(ts.frame, status)

        self.densities = {}
        for c in self.collectors:
            c.finish()
            self.densities[c.name] = c.Density()
        # should save precious files!!!
        return self.densities

    def DensityWithBulk(self, density_unit='water', solvent_threshold=numpy.e, bulk_threshold=0.6):
        """Return a solvent density with bulk site inserted.

        DensityWithBulk(self, solvent_threshold=2.72, bulk_threshold=0.6) --> Density

        Only works if two densities were generated that are named 'solvent' and
        'bulk' (this is the default for the *mode* = "all" keyword for
        :class:`DensityCreator`.)

        :Arguments:
          *density_unit*
              Measure density in multiples of this unit; possible values are
              'Molar', 'nm', 'Angstrom', or the density at standard conditions
              of 'water' (experimental value), 'TIP3P', 'TIP4P', 'SPC' ['water']
          *solvent_threshold*
              Hydration sites are considered as regions of density > this
              threshold; it is assumed to be given in the *density_unit*.
          *bulk_threshold*
              The bulk site is the largest region with a density >
              *bulk_threshold*; in order to avoid overlap with the hydration
              sites, it is necessary to use a special selection for the solvent
              that excludes it from the vicinity of the solute.

        .. SeeAlso:: This method uses meth:`hop.sitemap.Density.map_sites` and
                     meth:`hop.sitemap.Density.site_insert_bulk`.
        """
        if len(self.densities) != 2:
            errmsg = "DensityCreator.DensityWithBulk(): Need exactly two densities ('solvent' and 'bulk')."
            logger.fatal(errmsg)
            raise MissingDataError(errmsg)
        try:
            solvent = self.densities['solvent']
            bulk = self.densities['bulk']
        except KeyError:
            errmsg = "Need a 'solvent' and a 'bulk' density in %s.densities" % self.__class__.__name__
            logger.fatal(errmsg)
            raise MissingDataError(errmsg)
        logger.debug("DensityWithBulk: solvent_threshold = %r", solvent_threshold)
        logger.debug("DensityWithBulk: bulk_threshold = %r", bulk_threshold)
        solvent.convert_density(density_unit)
        solvent.map_sites(solvent_threshold)
        bulk.convert_density(density_unit)
        bulk.map_sites(bulk_threshold)

        # ye olde bulk-hack....
        solvent.site_insert_bulk(bulk)

        # should really save
        # solvent.save()

        return solvent



def density_from_dcd(*args,**kwargs):
    import warnings
    warnings.warn("density_from_dcd() is deprecated and will be removed. "
                  "Use density_from_trajectory().", category=DeprecationWarning)
    return density_from_trajectory(*args,**kwargs)

def density_from_trajectory(*args,**kwargs):
    """Create a density grid from a trajectory.

       density_from_trajectory(PSF, DCD, delta=1.0, atomselection='name OH2', ...) --> density

    or

       density_from_trajectory(PDB, XTC, delta=1.0, atomselection='name OH2', ...) --> density

    :Arguments:
      psf/pdb/gro
            topology file
      dcd/xtc/trr/pdb
            trajectory; if reading a single PDB file it is sufficient to just provide it
            once as a single argument

    :Keywords:
      atomselection
            selection string (MDAnalysis syntax) for the species to be analyzed
            ["name OH2"]
      delta
            approximate bin size for the density grid in Angstroem (same in x,y,z)
            (It is slightly adjusted when the box length is not an integer multiple
            of delta.) [1.0]
      metadata
            dictionary of additional data to be saved with the object
      padding
            increase histogram dimensions by padding (on top of initial box size)
            in Angstroem [2.0]
      soluteselection
            MDAnalysis selection for the solute, e.g. "protein" [``None``]
      cutoff
            With *cutoff*, select '<atomsel> NOT WITHIN <cutoff> OF <soluteselection>'
            (Special routines that are faster than the standard AROUND selection) [0]
      verbosity: int
            level of chattiness; 0 is silent, 3 is verbose [3]

    :Returns: :class:`hop.sitemap.Density`

    :TODO:
      * Should be able to also set skip and start/stop for data collection.

    .. Note::
        * In order to calculate the bulk density, use

              atomselection='name OH2',soluteselection='protein and not name H*',cutoff=3.5

          This will select water oxygens not within 3.5 A of the protein heavy atoms.
          Alternatively, use the VMD-based  :func:`density_from_volmap` function.
        * The histogramming grid is determined by the initial frames min and max.
        * metadata will be populated with psf, dcd, and a few other items.
          This allows more compact downstream processing.

    .. SeeAlso:: docs for :func:`density_from_Universe` (defaults for kwargs are defined there).
    """
    import MDAnalysis
    return density_from_Universe(MDAnalysis.Universe(*args),**kwargs)

def density_from_Universe(universe,delta=1.0,atomselection='name OH2',
                          metadata=None,padding=2.0,cutoff=0,soluteselection=None,verbosity=3,
                          use_kdtree=True, **kwargs):
    """Create a density grid from a MDAnalysis.Universe object.

      density_from_dcd(universe, delta=1.0, atomselection='name OH2', ...) --> density

    :Arguments:
      universe
            :class:`MDAnalysis.Universe` object with a trajectory

    :Keywords:
      atomselection
            selection string (MDAnalysis syntax) for the species to be analyzed
            ["name OH2"]
      delta
            approximate bin size for the density grid in Angstroem (same in x,y,z)
            (It is slightly adjusted when the box length is not an integer multiple
            of delta.) [1.0]
      metadata
            dictionary of additional data to be saved with the object
      padding
            increase histogram dimensions by padding (on top of initial box size)
            in Angstroem [2.0]
      soluteselection
            MDAnalysis selection for the solute, e.g. "protein" [``None``]
      cutoff
            With *cutoff*, select '<atomsel> NOT WITHIN <cutoff> OF <soluteselection>'
            (Special routines that are faster than the standard AROUND selection) [0]
      verbosity: int
            level of chattiness; 0 is silent, 3 is verbose [3]
      parameters
            dict with some special parameters for :class:`~hop.sitemap.Density` (see doc)
      kwargs
            metadata, parameters are modified and passed on to :class:`~hop.sitemap.Density`

    :Returns: :class:`hop.sitemap.Density`

    """
    try:
        universe.selectAtoms('all')
        universe.trajectory.ts
    except AttributeError:
        errmsg = "density_from_Universe(): The universe must be a proper MDAnalysis.Universe instance."
        logger.fatal(errmsg)
        raise TypeError(errmsg)
    u = universe
    if cutoff > 0 and soluteselection is not None:
        # special fast selection for '<atomsel> not within <cutoff> of <solutesel>'
        notwithin_coordinates = notwithin_coordinates_factory(u,atomselection,soluteselection,cutoff,use_kdtree=use_kdtree)
        def current_coordinates():
            return notwithin_coordinates()
    else:
        group = u.selectAtoms(atomselection)
        def current_coordinates():
            return group.coordinates()

    coord = current_coordinates()
    logger.info("Selected %d atoms out of %d atoms (%s) from %d total.",
                coord.shape[0],len(u.selectAtoms(atomselection)),atomselection,len(u.atoms))

    # mild warning; typically this is run on RMS-fitted trajectories and
    # so the box information is rather meaningless
    box,angles = u.trajectory.ts.dimensions[:3], u.trajectory.ts.dimensions[3:]
    if tuple(angles) <> (90.,90.,90.):
        logger.warn("Non-orthorhombic unit-cell --- make sure that it has been remapped properly!")

    # Make the box bigger to avoid as much as possible 'outlier'. This
    # is important if the sites are defined at a high density: in this
    # case the bulk regions don't have to be close to 1 * n0 but can
    # be less. It's much more difficult to deal with outliers.  The
    # ideal solution would use images: implement 'looking across the
    # periodic boundaries' but that gets complicate when the box
    # rotates due to RMS fitting.
    smin = numpy.min(coord,axis=0) - padding
    smax = numpy.max(coord,axis=0) + padding

    BINS = fixedwidth_bins(delta, smin, smax)
    arange = zip(BINS['min'],BINS['max'])
    bins = BINS['Nbins']

    # create empty grid with the right dimensions (and get the edges)
    grid,edges = numpy.histogramdd(numpy.zeros((1,3)),bins=bins,range=arange,normed=False)
    grid *= 0.0
    h = grid.copy()

    pm = CustomProgressMeter(u.trajectory.numframes, interval=10,
                             format="Histograming %(other)d atoms in frame %(step)5d/%(numsteps)d  [%(percentage)5.1f%%]\r")
    for ts in u.trajectory:
        coord = current_coordinates()
        if len(coord) == 0: continue
        h[:],edges[:] = numpy.histogramdd(coord, bins=bins, range=arange, normed=False)
        grid += h  # accumulate average histogram
        pm.echo(ts.frame, len(coord))
    numframes = u.trajectory.numframes / u.trajectory.skip
    grid /= float(numframes)

    # pick from kwargs
    metadata = kwargs.pop('metadata',{})
    metadata['psf'] = u.filename              # named psf for historical reasons
    metadata['dcd'] = u.trajectory.filename   # named dcd for historical reasons
    metadata['atomselection'] = atomselection
    metadata['numframes'] = numframes
    metadata['totaltime'] = round(u.trajectory.numframes * u.trajectory.delta * u.trajectory.skip_timestep \
                                  * constants.get_conversion_factor('time','AKMA','ps'), 3)
    metadata['dt'] = u.trajectory.delta * u.trajectory.skip_timestep * \
                     constants.get_conversion_factor('time','AKMA','ps')
    metadata['time_unit'] = 'ps'
    metadata['dcd_skip'] = u.trajectory.skip_timestep  # frames
    metadata['dcd_delta'] = u.trajectory.delta         # in AKMA
    if cutoff > 0 and soluteselection is not None:
        metadata['soluteselection'] = soluteselection
        metadata['cutoff'] = cutoff             # in Angstrom

    parameters = kwargs.pop('parameters',{})
    parameters['isDensity'] = False             # must override

    # all other kwargs are discarded

    # Density automatically converts histogram to density for isDensity=False
    g = Density(grid=grid,edges=edges,unit=dict(length='Angstrom'),
                parameters=parameters,metadata=metadata)
    logger.info("Histogram completed (initial density in Angstrom**-3)")

    return g


def notwithin_coordinates_factory(universe,sel1, sel2, cutoff, not_within=True, use_kdtree=True):
    """Generate optimized selection for '*sel1* not within *cutoff* of *sel2*'

    Example usage::
      notwithin_coordinates = notwithin_coordinates_factory(universe, 'name OH2','protein and not name H*',3.5)
      ...
      coord = notwithin_coordinates()        # changes with time step
      coord = notwithin_coordinates(cutoff2) # can use different cut off

    :Keywords:
      not_within
         True: selection behaves as 'not within' (As described above)
         False: selection is a <sel1> WITHIN <cutoff> OF <sel2>'
      use_kdtree
         True: use fast kd-tree based selections (requires new MDAnalysis >= 0.6)
         False: use distance matrix approach

    .. Note::
        * Periodic boundary conditions are NOT taken into account: the naive
          minimum image convention employed in the distance check is currently
          not being applied to remap the coordinates themselves, and hence it
          would lead to counts in the wrong region.
        * The selections are static and do not change with time steps.
    """
    # Benchmark of FABP system (solvent 3400 OH2, protein 2100 atoms) on G4 powerbook, 500 frames
    #                    cpu/s    relative   speedup       use_kdtree
    # distance matrix    633        1          1           False
    # AROUND + kdtree    420        0.66       1.5         n/a ('name OH2 around 4 protein')
    # manual + kdtree    182        0.29       3.5         True
    solvent = universe.selectAtoms(sel1)
    protein = universe.selectAtoms(sel2)
    if use_kdtree:
        # using faster hand-coded 'not within' selection with kd-tree
        import MDAnalysis.KDTree.NeighborSearch as NS
        import MDAnalysis.core.AtomGroup
        set_solvent = set(solvent)     # need sets to do bulk = allsolvent - selection
        if not_within is True:  # default
            def notwithin_coordinates(cutoff=cutoff):
                # must update every time step
                ns_w = NS.AtomNeighborSearch(solvent)  # build kd-tree on solvent (N_w > N_protein)
                solvation_shell = ns_w.search_list(protein,cutoff)  # solvent within CUTOFF of protein
                group = MDAnalysis.core.AtomGroup.AtomGroup(set_solvent - set(solvation_shell)) # bulk
                return group.coordinates()
        else:
            def notwithin_coordinates(cutoff=cutoff):
                # acts as '<solvent> WITHIN <cutoff> OF <protein>'
                # must update every time step
                ns_w = NS.AtomNeighborSearch(solvent)  # build kd-tree on solvent (N_w > N_protein)
                group = ns_w.search_list(protein,cutoff)  # solvent within CUTOFF of protein
                return group.coordinates()
    else:
        # slower distance matrix based (calculate all with all distances first)
        import MDAnalysis.core.distances
        dist = numpy.zeros((len(solvent),len(protein)),dtype=numpy.float64)
        box = None  # as long as s_coor is not minimum-image remapped
        if not_within is True:   # default
            compare = numpy.greater
            aggregatefunc = numpy.all
        else:
            compare = numpy.less_equal
            aggregatefunc = numpy.any
        def notwithin_coordinates(cutoff=cutoff):
            s_coor = solvent.coordinates()
            p_coor = protein.coordinates()
            # Does water i satisfy d[i,j] > r for ALL j?
            d = MDAnalysis.distances.distance_array(s_coor,p_coor,box=box,result=dist)
            return s_coor[aggregatefunc(compare(d,cutoff), axis=1)]
    return notwithin_coordinates



def density_from_volmap(psf,dcd,dx=None,delta=1.0,atomselection='name OH2',
                        metadata=None,verbosity=3,**kwargs):
    """Create the density using VMD's VolMap plugin and an intermediate dx file.

    density_from_volmap(psf,dcd,dx,delta=1.0,atomselection='name OH2',
          metadata=dict(psf=psf,dcd=dcd,system='I-FABP apo')

    psf     Charmm psf topology file
    dcd     Charmm trajectory
    dx      dx file that holds the density generated by VolMap; if not provided a
            temporary file is created (and deleted)
    atomselection
            selection string (MDAnalysis syntax) for the species to be analyzed
    delta   approximate bin size for the density grid (same in x,y,z)
            (It is slightly adjusted when the box length is not an integer multiple
            of delta.)
    verbosity=int  level of chattiness; 0 is silent, 3 is verbose

    **kwargs   metadata and parameters are passed to Density(), everything else to VolMap

     Density args:
     metadata       dictionary of additional data to be saved with the object
     parameters     dict of special Density parameters (see Density() doc)

     VolMap args:
     load_new       True: load psf and dcd into VMD. False: use psf and dcd already
                   loaded into VMD (default is True)

    Returns a Density object.

    """
    import hop.external
    import tempfile

    set_verbosity(verbosity)  # set to 0 for no messages

    # pick Density parameters from kwargs
    parameters = kwargs.pop('parameters',{})
    metadata = kwargs.pop('metadata',{})
    # ... everything else is VolMap args

    remove_dx = False
    if not dx:
        remove_dx = True
        dxhandle,dx = tempfile.mkstemp('.dx')
        logger.debug("Using tempororary dx file %r.", dx)

    logger.info("Connecting to VMD (ignore 'error: uncaptured python exception')")
    vmd = hop.external.VMD()
    logger.info("VolMap calculates the density. This takes a while...")
    vmd.volmap(psf,dcd,dx,delta=delta,atomselection=atomselection,**kwargs)

    metadata['psf'] = psf
    metadata['dcd'] = dcd
    metadata['vmd_dx'] = dx
    metadata['atomselection'] = atomselection

    parameters['isDensity'] = True             # must override

    logger.info("Building density object from dx file '%(dx)s'...\n", vars())
    g = Density(dxfile=dx,unit=dict(length='Angstrom',density='Angstrom'),
                parameters=parameters,metadata=metadata)

    if remove_dx:
        os.remove(dx)
    return g

def Bfactor2RMSF(B):
    """Atomic root mean square fluctuation (in Angstrom) from the crystallographic B-factor

    B = [(8*PI**2)/3] * (RMSF)**2

    Willis & Pryor, Thermal vibrations in crystallography, Cambridge
    Univ. Press, 1975
    """
    return numpy.sqrt(3.*B/8.)/numpy.pi

class PDBDensity(Density):
    __doc__ = """Density with additional information about original crystal structure.

This is simply the Density class (see below) enhanced by the add_xray2psf(),
W(), and Wequiv() methods.

Note that later analysis often ignores the site with the bulknumber by default
so one should (after computing a site map) also insert an empty bulk site:

  # canonical way to build a PDBDensity
  # (builds the sitepa at threshold and inserts a pseudo bulk site)
  xray = BfactorDensityCreator(...).PDBDensity(threshold)

  # rebuild site map
  xray.map_sites(threshold) # map sites at density cutoff threshold
  xray.site_insert_nobulk() # insert 'fake' bulk site at position SITELABEL['bulk']

  # find X-ray waters that correspond to a site in another density Y:
  # (1) build the list of equivalence sites, using the x-ray density as reference
  Y.find_equivalence_sites(xray)     # also updates equiv-sites in xray!
  # (2) look at the matches in xray
  xray.Wequiv()        TODO: not working yet


""" + 60*"-" + "\nDensity Class\n\n" + Density.__doc__

    # will probably break under multiple inheritance but I haven't figured out how to use super here
    _saved_attributes = Density._saved_attributes + ['_xray2psf', '_psf2xray']

    def add_xray2psf(self,pdbfile,regex=r'\s*W\s*|HOH|WAT|.*TIP.*|.*SPC.*'):
        """Add translation table between sequential psf numbering and original pdb numbering for water.

        D.add_xray2psf(pdbfilename)

        The original pdb is read and all water molecules are sequentially mapped
        to the water molecules in the psf (without any checks). The pdb is read
        and analyzed using Bio.PDB.

        pdbfilename    Original crystallographic pdb file
        regex          extended regular expression to detect water   residues

        """
        import re
        import Bio.PDB
        water = re.compile(regex)
        parser = Bio.PDB.PDBParser()
        m = parser.get_structure('0UNK',pdbfile)
        s = m[0]
        # number waters sequentially and store the pdb resid
        self._xray2psf = dict([(resid_xray,resid_psf+1) for resid_psf,resid_xray in
                               enumerate([r.id[1] for r in s.get_residues() if water.match(r.resname)])
                               ])
        self._psf2xray = dict([(resid_psf,resid_xray) for resid_xray,resid_psf in self._xray2psf.items()])
    def _check_site_resid_match(self):
        return len(self._xray2psf) == len(self.site_labels('sites',exclude='equivalencesites'))

    def W(self,N,returntype="auto",format=False):
        """Returns the resid of water N.

        If returntype == 'psf' then N is interpreted as the resid in the
        x-ray crystal structure (or original pdb file) and a resid N' in the
        psf is returned.

        If returntype == 'xray' then N is a resid in the psf and the
        corresponding crystal structure water is returned. This is
        useful to label water molecules by their published identifier,
        eg 'W128'.

        If the returntype is set to 'auto' and N starts with a W (eg
        'W128') then it is assumed to be a crystal water and the
        returntype is automatically set to psf, otherwise it acts like
        'xray'.

        :Arguments:
        N              resid of molecule (can be an iterable)
        returntype     'auto' | 'psf' | 'xray'
        format         False: return a integer number
                       True: default string (either "WN'" for x-ray or "#N'" for psf)
                       python format string: if the string contains %(resid)d then the string
                                             will be used as a format, otherwise the bare number
                                             is returned without raising an error
        """
        if returntype not in ("auto","psf","xray"):
            raise ValueError("returntype must be one of 'psf' or 'xray'")
        result = numpy.array([self._getN(_N,returntype=returntype,format=format) for _N in asiterable(N)])
        if not iterable(N):
            return result[0]
        else:
            return result

    def Wequiv(self,format=True):
        """Return a list of the PDB resids of the equivalent sites.

        array = Wequiv(format=True)

        format        True: array of identifiers 'Wnn'
                      False: array of integers
                      string: python format string; %(resid)d is replaced
        """
        return self.W(self.site_labels('subsites'),format=format)

    def _getN(self,N,returntype='auto',format=False):
        if returntype is 'auto':
            _Nstring = str(N).upper()
            returntype = "xray"
            if _Nstring.startswith('W'):
                # automagically do the right thing
                returntype = "psf"
                N = int(_Nstring[1:])
            elif _Nstring.startswith('#'):
                N = int(_Nstring[1:])
        if returntype == "psf":
            return self._Wpsf(N,format=format)
        elif returntype == "xray":
            return self._Wxray(N,format=format)

    def _Wpsf(self,resid_xray,format=False):
        """Returns resid in psf of crystallographic water W(resid_xray)."""
        try:
            resid = self._xray2psf[resid_xray]
        except KeyError:
            raise ValueError("No residue number %(resid_xray)d in x-ray structure." % vars())
        except AttributeError:
            raise MissingDataError("Add the xray -> psf translation table with add_xray2psf() first.")
        return self._Wformatter(resid,format=format,typechar='#')

    def _Wxray(self,resid_psf,format=False):
        """Returns the crystal structure resid of water resid_psf in the psf."""
        try:
            resid = self._psf2xray[resid_psf]
        except KeyError:
            raise ValueError("No residue number %(resid_psf)d in psf." % vars())
        except AttributeError:
            raise MissingDataError("Add the psf -> x-ray translation table with add_xray2psf() first.")
        return self._Wformatter(resid,format=format,typechar='W')

    def _Wformatter(self,resid,format=False,typechar=None):
        # no error checks, only call from wrappers
        if format is True and typechar is not None:
            default_format = str(typechar)+'%(resid)d'
            return  default_format % vars()
        elif str(format).find('%(resid)') >= 0:
            return str(format) % vars()
        else:
            return resid

    def site_insert_nobulk(self):
        """Insert an empty bulk site for cases when this is convenient."""
        class Nobulk:
            def __init__(self,dens):
                # copy the attributes that are checked in Density.site_insert_bulk()
                self.map = numpy.empty_like(dens.map)
                self.unit = dens.unit
                self.P = {'threshold': None}
                # minimum empty sites 'list'
                self.sites = {SITELABEL['bulk']: ()}  # normally a list but use a dict :-)
        self.site_insert_bulk(Nobulk(self))

    def equivalence_sites(self,format=True):
        """All equivalence sites (if defined) together with crystallographic water labels.

        recarray <-- equivalence_sites(self,format=True)

        The numpy.recarray has columns
            equivalence_label        the integer label of the equivalence site
            equivalence_name         the name, a string
            xray                     the identifier of the X-ray water

        equivalence_label and equivalence_name are identical between the densities from
        which the equivalence sites were computed. The xray identifier is specific for the
        structure; by default it is a string such as 'W135'.

        format        True: print 'W<N>' identifier
                      False: integer <N>
                      (see W() for more possibilities)

        BUG: THIS IS NOT WORKING AS THOUGHT BECAUSE THERE IS NO 1-1
        MAPPING BETWEEN WATER MOLECULES AND SITES AND BECAUSE SITES
        ARE NOT NUMBERED IN THE SAME ORDER AS THE WATER MOLECULES

        TODO: The proper way to do this is to find all water molecules
        within a cutoff of each grid cell that belongs to a site and
        then store all the waters as the string name of the site.

        """

        records = []
        equiv = self.subsites_of(self.site_labels('equivalencesites'))
        for equivlabel,subsites in equiv.items():
            records.append((equivlabel, self._equivlabel2equivname(equivlabel)[0],
                            self.W(self.site2resid(subsites.label[0]), returntype='xray', format=format)))
        return numpy.rec.fromrecords(records,names="equivalence_label,equivalence_name,xray")

    def site2resid(self,sitelabel):
        """Returns the resid of the particle that provided the density  for the site.
        """
        raise NotImplementedError('site2resid mapping is not working yet')

def print_combined_equivalence_sites(target,reference):
    """Tabulate equivalence sites of target against the reference.

    BUG: THIS IS NOT WORKING (because the assignment sites <--> waters
    is broken)
    """
    raise NotImplementedError('THIS IS NOT WORKING (because the assignment sites <--> waters is broken')

    eqs_r = reference.equivalence_sites()
    eqs_t = target.equivalence_sites()
    eqs_r.equivalence_name.sort()
    sorted_t = eqs_t[eqs_t.equivalence_label == eqs_r.equivalence_label]
    _header =  "%3s %4s   %-5s %-6s" % ('i','name', 'ref', 'target')
    def ___(): print '-' * len(_header)
    ___()
    print _header
    ___()
    for (l,n,x1),(l,n,x2) in zip(eqs_r,sorted_t):
        print "%3d %4s   %-5s  %-5s" % (l,n,x1,x2)
    ___()


class BfactorDensityCreator(object):
    """Create a density grid from a pdb file using MDAnalysis.

      dens = BfactorDensityCreator(psf,pdb,...).PDBDensity()

    The main purpose of this function is to convert crystal waters in
    an X-ray structure into a density so that one can compare the
    experimental density with the one from molecular dynamics
    trajectories. Because a pdb is a single snapshot, the density is
    estimated by placing Gaussians of width sigma at the position of
    all selected atoms.

    Sigma can be fixed or taken from the B-factor field, in which case
    sigma is taken as sqrt(3.*B/8.)/pi.

    TODO:

    * Make Gaussian convolution more efficient (at least for same
      sigma) because right now it is VERY slow (which may be
      acceptable if one only runs this once)
    * Using a temporary Creator class with the PDBDensity() helper
      method is clumsy (but was chosen as to keep the PDBDensity class
      clean and __init__ compatible with Density).
    """
    def __init__(self, psf,pdb,delta=1.0,atomselection='name OH2',
                metadata=None,padding=4.0, sigma=None, verbosity=3):
        """Construct the density from psf and pdb and the atomselection.

        DC = BfactorDensityCreator(psf, pdb, delta=<delta>, atomselection=<MDAnalysis selection>,
                                  metadata=<dict>, padding=2, sigma=None)
        density = DC.PDBDensity()

        psf     Charmm psf topology file
        pdb     PDB file
        atomselection
                selection string (MDAnalysis syntax) for the species to be analyzed
        delta   approximate bin size for the density grid (same in x,y,z)
                (It is slightly adjusted when the box length is not an integer multiple
                of delta.)
        metadata
                dictionary of additional data to be saved with the object
        padding increase histogram dimensions by padding (on top of initial box size)
        sigma   width (in Angstrom) of the gaussians that are used to build up the
                density; if None then uses B-factors from pdb
        verbosity=int  level of chattiness; 0 is silent, 3 is verbose

        For assigning X-ray waters to MD densities one might have to use a sigma
        of about 0.5 A to obtain a well-defined and resolved x-ray water density
        that can be easily matched to a broader density distribution.

        """
        from MDAnalysis import Universe
        set_verbosity(verbosity)  # set to 0 for no messages
        u = Universe(psf,pdbfilename=pdb)
        group = u.selectAtoms(atomselection)
        coord = group.coordinates()
        logger.info("BfactorDensityCreator: Selected %d atoms (%s) out of %d total.",
                    coord.shape[0],atomselection,len(u.atoms))
        smin = numpy.min(coord,axis=0) - padding
        smax = numpy.max(coord,axis=0) + padding

        BINS = fixedwidth_bins(delta, smin, smax)
        arange = zip(BINS['min'],BINS['max'])
        bins = BINS['Nbins']

        # get edges by doing a fake run
        grid,self.edges = numpy.histogramdd(numpy.zeros((1,3)),
                                            bins=bins,range=arange,normed=False)
        self.delta = numpy.diag(map(lambda e: (e[-1] - e[0])/(len(e)-1), self.edges))
        self.midpoints = map(lambda e: 0.5 * (e[:-1] + e[1:]), self.edges)
        self.origin = map(lambda m: m[0], self.midpoints)
        numframes = 1

        if sigma is None:
            # histogram individually, and smear out at the same time
            # with the appropriate B-factor
            if numpy.any(group.bfactors == 0.0):
                wmsg = "BfactorDensityCreator: Some B-factors are Zero."
                warnings.warn(wmsg, category=hop.MissingDataWarning)
                logger.warn(wmsg)
            rmsf = Bfactor2RMSF(group.bfactors)
            grid *= 0.0  # reset grid
            self.g = self._smear_rmsf(coord,grid,self.edges,rmsf)
        else:
            # histogram 'delta functions'
            grid,self.edges = numpy.histogramdd(coord,bins=bins,range=arange,normed=False)
            logger.info("Histogrammed %6d atoms from pdb.", len(group.atoms))
            # just a convolution of the density with a Gaussian
            self.g = self._smear_sigma(grid,sigma)

        try:
            metadata['psf'] = psf
        except TypeError:
            metadata = dict(psf=psf)
        metadata['pdb'] = pdb
        metadata['atomselection'] = atomselection
        metadata['numframes'] = numframes
        metadata['sigma'] = sigma
        self.metadata = metadata

        # Density automatically converts histogram to density for isDensity=False
        logger.info("BfactorDensityCreator: Histogram completed (initial density in Angstrom**-3)\n")


    def PDBDensity(self,threshold=None):
        """Returns a PDBDensity object.

        The PDBDensity is a Density with a xray2psf translation table;
        it has also got an empty bulk site inserted (so that any
        further analysis which assumes that site number 1 is the bulk)
        does not discard a valid site.

        threshold      Use the given threshold to generate the graph; the threshold
                       is assumed to be in the same units as the density.
                       None: choose defaults (1.0 if bfactors were used, 1.3 otherwise)

        """
        d = PDBDensity(grid=self.g,edges=self.edges,unit=dict(length='Angstrom'),
                       parameters=dict(isDensity=False),metadata=self.metadata)
        d.add_xray2psf(d.metadata['pdb'])     # pdb filename is recorded in metadata
        d.convert_density('water')
        if threshold is None:
            if self.metadata['sigma'] is None:
                threshold = 1.0
            else:
                threshold = 1.3
        d.map_sites(threshold)
        d.site_insert_nobulk()   # fake bulk site
        if not d._check_site_resid_match():
            wmsg = "BfactorDensityCreator.PDBDensity(): "\
                "There are different numbers of water molecules (%d) and sites (%d). " \
                "Site <-> water matching will not work."  % \
                (len(d._xray2psf), len(d.site_labels('sites', exclude='equivalencesites')))
            logger.warn(wmsg)
            warnings.warn(wmsg, category=InconsistentDataWarning)
        return d

    def _smear_sigma(self,grid,sigma):
        # smear out points
        # (not optimized -- just to test the principle; faster approach could use
        # convolution of the whole density with a single Gaussian via FFTs:
        # rho_smeared = F^-1[ F[g]*F[rho] ]
        g = numpy.zeros(grid.shape)   # holds the smeared out density
        pos = numpy.where(grid <> 0)  # position in histogram (as bin numbers)
        pm = CustomProgressMeter(len(pos[0]), interval=1, offset=1,
                                 format="Smearing out water position %(step)4d/%(numsteps)5d with RMSF %(other)4.2f A [%(percentage)5.2f%%\r")
        for iwat in xrange(len(pos[0])): # super-ugly loop
            p = tuple([wp[iwat] for wp in pos])
            g += grid[p] * \
                numpy.fromfunction(self._gaussian,grid.shape,dtype=numpy.int,
                                   p=p,sigma=sigma)
            pm.echo(iwat, sigma)
        return g

    def _smear_rmsf(self,coordinates,grid,edges,rmsf):
        # smear out each water with its individual Gaussian
        # (slower than smear_sigma)
        g = numpy.zeros(grid.shape)   # holds the smeared out density
        N,D = coordinates.shape
        pm = CustomProgressMeter(N, interval=1, offset=1,
                                 format="Smearing out water position %(step)4d/%(numsteps)5d with RMSF %(other)4.2f A [%(percentage)5.2f%%\r")
        for iwat,coord in enumerate(coordinates):
            g += numpy.fromfunction(self._gaussian_cartesian,grid.shape,dtype=numpy.int,
                                    c=coord,sigma=rmsf[iwat])
            pm.echo(iwat, rmsf[iwat])
        return g

    def _gaussian(self,i,j,k,p,sigma):
        # i,j,k can be numpy arrays
        # p is center of gaussian as grid index, sigma its width (in A)
        x = self.delta[0,0]*(i - p[0])  # in Angstrom
        y = self.delta[1,1]*(j - p[1])
        z = self.delta[2,2]*(k - p[2])
        return (2*numpy.pi*sigma)**(-1.5) * numpy.exp(-(x*x+y*y+z*z)/(2*sigma*sigma))

    def _gaussian_cartesian(self,i,j,k,c,sigma):
        # i,j,k can be numpy arrays
        # c is center of gaussian in cartesian coord (A), sigma its width (in A)
        x = self.origin[0] + self.delta[0,0]*i - c[0]  # in Angstrom
        y = self.origin[1] + self.delta[1,1]*j - c[1]
        z = self.origin[2] + self.delta[2,2]*k - c[2]
        return (2*numpy.pi*sigma)**(-1.5) * numpy.exp(-(x*x+y*y+z*z)/(2*sigma*sigma))


