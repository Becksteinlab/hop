# special hop classes that have not been moved into MDAnalysis.analysis.density
# Already changed msg --> logger and various other bits and pieces.

import numpy  # need v >= 1.0
import sys
import os,os.path,errno
import cPickle
import warnings

from gridData import Grid, OpenDX    # http://github.com/orbeckst/GridDataFormats 

import MDAnalysis
from MDAnalysis.core.util import fixedwidth_bins, iterable, asiterable
from MDAnalysis import NoDataError

import logging
logger = logging.getLogger("MDAnalysis.analysis.density")


class DensityCollector(object):
    """Collect subsequent coordinate frames to build up a :class:`Density`."""
    use_kdtree = True
    def __init__(self, name, universe, **kwargs):
        self.name = name
        try:
            universe.selectAtoms('all')
            universe.trajectory.ts
        except AttributeError:
            raise TypeError("The universe must be a proper MDAnalysis.Universe instance.")        
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
        logger.info("%-10s: Selected %d atoms out of %d atoms (%s) from %d total." %
                    (self.name, coord.shape[0],len(u.selectAtoms(self.atomselection)),
                     self.atomselection,len(u.atoms)))

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
        """Return a :class:`Density` from the data."""
        if not hasattr(self, 'grid'):
            raise NoDataError("No data for density available. Run collect() first.")
        u = self.universe
        metadata = self.metadata
        metadata['collector'] = self.name
        metadata['collector_mode'] = self.mode
        metadata['topology'] = u.filename
        metadata['trajectory'] = u.trajectory.filename
        metadata['atomselection'] = self.atomselection
        metadata['numframes'] = u.trajectory.numframes
        metadata['dt'] = u.trajectory.dt    # in ps for default MDAnalysis
        # totaltime should be in MDAnalysis!
        metadata['totaltime'] = round(u.trajectory.numframes * metadata['dt'] * u.trajectory.skip_timestep, 3)
        metadata['time_unit'] = MDAnalysis.core.flags['time_unit']  # just to make sure we know it...
        metadata['skip_timestep'] = u.trajectory.skip_timestep  # frames
        metadata['delta'] = u.trajectory.delta         # in native units (?)
        if self.mode == 'BULK':
            metadata['soluteselection'] = self.soluteselection
            metadata['cutoff'] = self.cutoff             # in Angstrom

        parameters = self.parameters
        parameters['isDensity'] = False             # must override

        # Density automatically converts histogram to density for isDensity=False
        g = Density(grid=self.grid, edges=self.edges,
                    unit=dict(length=MDAnalysis.core.flags['length_unit']),
                    parameters=parameters, metadata=metadata)    
        logger.info("%-10s: Histogram completed (initial density in %s**-3)" % (self.name, MDAnalysis.core.flags['length_unit'])) 
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
            raise ValueError("mode must be one of %r, not %r" % (self.modes, self.mode))
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
        u = self.universe
        for ts in u.trajectory:
            status = ["Histograming"]
            for c in self.collectors:
                natoms = c.collect()
                status.append("%s=%d" % (c.name, natoms))

            if u.trajectory.ts.frame % 10 == 0 or \
                    u.trajectory.ts.frame == u.trajectory.numframes:
                message = " ".join(status)
                message += " atoms in frame %5d/%d  [%5.1f%%]\r" % (
                    u.trajectory.ts.frame, 
                    u.trajectory.numframes,
                    100.0*u.trajectory.ts.frame/u.trajectory.numframes)
                print message
        print
        self.densities = {}
        for c in self.collectors:
            c.finish()
            self.densities[c.name] = c.Density()
        # should save precious files!!!
        return self.densities

    def DensityWithBulk(self, density_unit='water', solvent_threshold=2.72, bulk_threshold=0.6):
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
            raise NoDataError("Need exactly two densities.")
        try:
            solvent = self.densities['solvent']
            bulk = self.densities['bulk']
        except KeyError:
            raise NoDataError("Need a 'solvent' and a 'bulk' density in %s.densities" %
                                       self.__class__.__name__)
        solvent.convert_density(density_unit)
        solvent.map_sites(solvent_threshold)
        bulk.convert_density(density_unit)
        bulk.map_sites(bulk_threshold)

        # ye olde bulk-hack....
        solvent.site_insert_bulk(bulk)

        # should really save
        # solvent.save()

        return solvent
