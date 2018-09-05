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
import MDAnalysis.analysis.density
from MDAnalysis.analysis.density import notwithin_coordinates_factory


from . import constants
from .exceptions import MissingDataError, InconsistentDataWarning
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
            universe.select_atoms('all')
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
            group = u.select_atoms(self.atomselection)
            self.current_coordinates = lambda : group.positions
            self.mode = "SOLVENT"
        coord = self.current_coordinates()
        logger.info("%-10s: Selected %d atoms out of %d atoms (%s) from %d total.",
                    self.name, coord.shape[0],len(u.select_atoms(self.atomselection)),
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

    def finish(self, n_frames):
        if self.isComplete():
            return
        self.grid /= float(n_frames)
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
        metadata['n_frames'] = u.trajectory.n_frames
        metadata['dt'] = u.trajectory.dt    # in ps for default MDAnalysis
        metadata['totaltime'] = u.trajectory.totaltime
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
        self.universe = MDAnalysis.as_Universe(*args, **universe_kwargs)
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
        pm = CustomProgressMeter(u.trajectory.n_frames, interval=10,
                                 format="Histogramming %(other)s atoms in frame %(step)5d/%(numsteps)d  [%(percentage)5.1f%%]\r")

        for ts in u.trajectory:
            status = ""
            for c in self.collectors:
                natoms = c.collect()
                status += ("%s=%d " % (c.name, natoms))
            pm.echo(ts.frame, status)

        self.densities = {}
        for c in self.collectors:
            c.finish(u.trajectory.n_frames)  # adjust if we implement trajectory slicing
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
              'Molar', 'nm^{-3}', 'Angstrom^{-3}', or the density at standard conditions
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

def density_from_Universe(*args, **kwargs):
    """Create a :class:`hop.sitemap.Density from a :class:`Universe`.

    .. SeeAlso::
       :func:`MDAnalysis.analysis.density.density_from_Universe` for
       all parameters and :func:`density_from_trajectory` for a
       convenience wrapper.

    """
    D = MDAnalysis.analysis.density.density_from_Universe(*args, **kwargs)
    return Density(grid=D.grid, edges=D.edges, units=D.units,
                   parameters=D.parameters, metadata=D.metadata)


def density_from_trajectory(*args, **kwargs):
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

    .. SeeAlso:: docs for
                 :func:`MDAnalysis.analysis.density.density_from_Universe`
                 (defaults for kwargs are defined there).

    """
    return density_from_Universe(MDAnalysis.Universe(*args),**kwargs)


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
        class Nobulk():
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


class BfactorDensityCreator(MDAnalysis.analysis.density.BfactorDensityCreator):
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

    .. SeeAlso::
       * :mod:`MDAnalysis.analysis.density`
       * :class:`PDBDensity`

    """

    def PDBDensity(self, threshold=None):
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



