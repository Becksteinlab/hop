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
Defining solvation sites --- :mod:`hop.sitemap`
===============================================

Histogram positions of particles from a MD trajectory on a
grid. Calculate the density, change units (both of the grid and of the
density), save the density, export into 3D visualization formats,
manipulate the density as a numpy array.
"""
from __future__ import absolute_import

import sys
import os, os.path
import errno
import cPickle
import pickle
import warnings

import numpy                   # need v >= 1.0
import networkx as NX          # https://networkx.lanl.gov/
from gridData import OpenDX    # https://github.com/MDAnalysis/GridDataFormats

from . import constants
from .constants import SITELABEL
from .exceptions import (MissingDataError, MissingDataWarning,
                         InconsistentDataWarning, OverwriteWarning)
from . import utilities
from .utilities import msg,set_verbosity,get_verbosity, flatten, sorted, \
     DefaultDict, fixedwidth_bins, iterable, asiterable
from itertools import izip

class Grid(utilities.Saveable):
    """Class to manage a multidimensional grid object.

    The grid (Grid.grid) can be manipulated as a standard numpy
    array. Changes can be saved to a file using the save() method. The
    grid can be restored using the load() method or by supplying the
    filename to the constructor.

    The attribute Grid.metadata holds a user-defined dictionary that
    can be used to annotate the data. It is saved with save().

    The export(format='dx') method always exports a 3D object, the
    rest should work for an array of any dimension.
    """

    _saved_attributes = ['grid','edges','P','unit','metadata']

    parameters_default = {'isDensity':False}
    unit_default = {'length':'Angstrom','density':None}

    def __init__(self,grid=None,edges=None,filename=None,dxfile=None,
                 parameters=None,unit=None,metadata=None):
        """
        Create a Grid object from data.

        From a numpy.histogramdd():
          g = Grid(grid,edges)
        From files (created with Grid.save(<filename>):
          g = Grid(filename=<filename>)
        From a dx file:
          g = Grid(dxfile=<dxfile>)

        Arguments:

        grid       histogram or density and ...
        edges      list of arrays, the lower and upper bin edges along the axes
                   (both are output by numpy.histogramdd())
        filename   file name of a pickled Grid instance (created with
                   Grid.save(filename))
        dxfile     OpenDX file
        parameters dictionary of class parameters; saved with save()
                   isDensity  False: grid is a histogram with counts,
                              True: a density.
                              Applying Grid.make_density() sets it to True.
        unit       dict(length='Angstrom', density=None)
                   length:  physical unit of grid edges (Angstrom or nm)
                   density: unit of the density if isDensity == True or None
        metadata   a user defined dictionary of arbitrary values
                   associated with the density; the class does not touch
                   metadata[] but stores it with save()

        Returns:
        g          a Grid object

        If the input histogram consists of counts per cell then the
        make_density() method converts the grid to a physical
        density. For a probability density, divide it by grid.sum() or
        use normed=True right away in histogramdd().

        If grid, edges, AND filename are given then the
        extension-stripped filename is stored as the default filename.

        NOTE:

        * It is suggested to construct the Grid object from a
          histogram, to supply the appropriate length unit, and to use
          make_density() to obtain a density. This ensures that the
          length- and the density unit correspond to each other.

        TODO:
        * arg list is still messy
        * probability density not supported as a unit
        """

        parameters = DefaultDict(self.parameters_default,parameters)
        unit = DefaultDict(self.unit_default,unit)
        metadata = DefaultDict({},metadata)

        # First set attributes that may be overriden by reading from a file
        # using load().
        self._dxfile = dxfile
        self.P = parameters    # isDensity: set by make_density()
        self.metadata = metadata                 # use this to record arbitrary data
        self.unit = unit
        self._check_set_unit(unit)   # unit must be dict --- check here?

        if not (grid is None or edges is None):
            self.grid = numpy.asarray(grid)
            self.edges = edges
        elif filename:
            super(Grid,self).__init__(filename=filename)  # uses self.load(filename)
            self._update()          # crucial because not all attributes are pickled!
            return                  # get out right away to keep it clean
        elif dxfile:
            self.importdx(dxfile)   # calls __init__() again with data
            return                  # get out right away to keep it clean
        else:
            raise ValueError("no input data---eg (grid,edges) or filename---for grid given")
        self._update()

    def _update(self):
        """compute/update all derived data

        Can be called without harm (idem-potent); needed separately
        when the units are changed.

        origin  is the center of the cell with index 0,0,0
        """
        self.delta = numpy.diag(
            map(lambda e: (e[-1] - e[0])/(len(e)-1), self.edges) )
        self.midpoints = map(lambda e: 0.5 * (e[:-1] + e[1:]), self.edges)
        self.origin = map(lambda m: m[0], self.midpoints)

        # sanity checks
        if self.P['isDensity']:
            if not self.unit['density']:
                raise ValueError("For a density, unit['density'] must be set.")

    def _check_set_unit(self,u):
        """Check that all bindings {unit_type : value, ...} in the dict u are valid
        and set the object's unit attribute.
        """
        # all this unit crap should be a class...
        # Also see comments near definition of hop.constants.conversion_factor[].
        try:
            for unit_type,value in u.items():
                if value is None:   # check here, too iffy to use dictionary[None]=None
                    self.unit[unit_type] = None
                    continue
                try:
                    constants.conversion_factor[unit_type][value]
                    self.unit[unit_type] = value
                except KeyError:
                    raise ValueError('Unit '+str(value)+\
                                     ' of type '+str(unit_type)+' is not recognized.')
        except AttributeError:
            raise ValueError('"unit" must be a dictionary with keys "length" and "density".')
        # need at least length and density (can be None)
        if 'length' not in self.unit:
            raise ValueError('"unit" must contain a unit for "length".')
        if 'density' not in self.unit:
            self.unit['density'] = None

    def make_density(self):
        """Convert the grid (a histogram, counts in a cell) to a density (counts/volume).

        make_density()

        Note: (1) This changes the grid irrevocably.
              (2) For a probability density, manually divide by grid.sum().
        """
        # Make it a density by dividing by the volume of each grid cell
        # (from numpy.histogramdd, which is for general n-D grids)
        if self.P['isDensity']:
            raise RuntimeError("Grid is already a density.")

        dedges = map(numpy.diff,self.edges)
        D = len(self.edges)
        for i in xrange(D):
            shape = numpy.ones(D, int)
            shape[i] = len(dedges[i])
            self.grid /= dedges[i].reshape(shape)
        self.P['isDensity'] = True
        self.unit['density'] = self.unit['length']

    def convert_length(self,unit='Angstrom'):
        """Convert Grid object to the new unit:
        Grid.convert_length(<unit>)

        unit       Angstrom, nm

        This changes the edges but will not change the density; it is
        the user's responsibility to supply the appropriate unit if
        the Grid object is constructed from a density. It is suggested
        to start from a histogram and a length unit and use
        make_density().
        """
        if unit == self.unit['length']:
            return
        cvnfact = constants.get_conversion_factor('length',self.unit['length'],unit)
        self.edges = [x * cvnfact for x in self.edges]
        self.unit['length'] = unit
        self._update()        # needed to recalculate midpoints and origin

    def convert_density(self,unit='Angstrom'):
        """Convert the density to the physical units given by unit

        Grid.convert_to(<unit>)

        <unit> can be one of the following:

        Angstrom     particles/A**3
        nm           particles/nm**3
        SPC          density of SPC water at standard conditions
        TIP3P        ... see __water__['TIP3P']
        TIP4P        ... ...
        water        density of real water at standard conditions (0.997 g/cm**3)
        Molar        mol/l

        Note: (1) This only works if the initial length unit is provided.
              (2) Conversions always go back to unity so there can be rounding
                  and floating point artifacts for multiple conversions.

        There may be some undesirable cross-interactions with convert_length...
        """
        if not self.P['isDensity']:
            raise RuntimeError('The grid is not a density so convert_density() makes no sense.')
        if unit == self.unit['density']:
            return
        self.grid *= constants.get_conversion_factor('density',self.unit['density'],unit)
        self.unit['density'] = unit

    def centers(self):
        """Returns the coordinates of the centers of all grid cells as an iterator."""
        # crappy
        for idx in numpy.ndindex(self.grid.shape):
            # TODO: CHECK that this delta*(i,j,k) is really correct, even for non-diagonal delta
            # NOTE: origin is center of (0,0,0) (and already has index offset by 0.5)
            yield numpy.sum(self.delta * numpy.asarray(idx), axis=0) + self.origin


    def importdx(self,dxfile):
        """Initializes Grid from a OpenDX file."""

        dx = OpenDX.field(0)
        dx.read(dxfile)
        grid,edges = dx.histogramdd()
        self.__init__(grid=grid,edges=edges,parameters=self.P,unit=self.unit,
                      metadata=self.metadata,dxfile=dxfile)

    def export(self,filename=None,format="dx"):
        """export density to file using the given format; use 'dx' for visualization.

        export(filename=<filename>,format=<format>)

        The <filename> can be omitted if a default file name already
        exists for the object (e.g. if it was loaded from a file or it
        was saved before.) Do not supply the filename extension. The
        correct one will be added by the method.

        The default format for export() is 'dx'.

        Only implemented formats:

        dx        OpenDX (WRITE ONLY)
        python    pickle (use Grid.load(filename) to restore); Grid.save()
                  is simpler than export(format='python').
        """
        filename = self.filename(filename)
        if format == "dx":
            self._export_dx(filename)
        elif format == "python":
            self.save(filename)
        else:
            raise NotImplementedError("Exporting to format "+str(format)+\
                                      " is not implemented.")

    def _export_dx(self,filename):
        """Export the density grid to an OpenDX file. The file format
        is the simplest regular grid array and it is also understood
        by VMD's DX reader.

        For the file format see
        http://opendx.sdsc.edu/docs/html/pages/usrgu068.htm#HDREDF
        """

        filename = filename + '.dx'

        comments = [
            'OpenDX density file written by',
            '$Id$',
            'File format: http://opendx.sdsc.edu/docs/html/pages/usrgu068.htm#HDREDF',
            'Data are embedded in the header and tied to the grid positions.',
            'Unit of   length    is '+str(self.unit['length'])+'.',
            'Unit of   density   is '+str(self.unit['density'])+'-based.',
            'Data is written in C array order: In grid[x,y,z] the axis z is fastest',
            'varying, then y, then finally x, i.e. z is the innermost loop.']

        # write metadata in comments section
        if self.metadata:
            comments.append('Meta data stored with the python Grid object:')
        for k in self.metadata:
            comments.append('   '+str(k)+' = '+str(self.metadata[k]))
        comments.append('(Note: the VMD dx-reader chokes on comments below this line)')

        components = dict(
            positions = OpenDX.gridpositions(1,self.grid.shape,self.origin,self.delta),
            connections = OpenDX.gridconnections(2,self.grid.shape),
            data = OpenDX.array(3,self.grid),
            )
        dx = OpenDX.field('density',components=components,comments=comments)
        dx.write(filename)

    def __repr__(self):
        if self.P['isDensity']:
            grid_type = 'density'
        else:
            grid_type = 'histogram'
        try:
            fn = 'and default filename "'+str(self.filename())+'"'
        except ValueError:
            fn = ''
        return '<hop.sitemap.Grid '+grid_type+' with '+str(self.grid.shape)+' bins '+fn+'>'


# could we derive Density from MDAnalysis.analysis.density.Density?

class Density(Grid):
    """Class with an annotated density, i.e. additional information
    for each grid cell. Adds information about sites to the grid. A
    'site' consists of all connected grid cells with a density >=
    threshold.

    A site is defined as a set of at least 'MINsite' grid cells with
    density >= threshold that are located in each others' first and
    second nearest neighbour shell (of 26 cells, on the cubic
    lattice). A site is labelled by an integer 1..N. The interstitial
    is labelled '0'. By default, a site may consist of a single grid
    cell (MINsite == 1) but this can be changed by setting the
    parameter MINsite to another number >1.

    When neither grid nor edges are given then the density object can
    also be read from a pickled file (filename) or a OpenDX file
    (dxfile). In the latter case, care should be taken to properly set
    up the units and the isDensity parameter:

    >>> g = Density(dxfile='bulk.dx',parameters={'isDensity':True,'MINsite':1},
                    unit={'length':'Angstrom','density':'Angstrom'}, ....)

    Attributes:

        grid          density on a grid
        edges         the lower and upper edges of the grid cells along the
                      three dimensions of the grid
        map           grid with cells labeled as sites (after label_sites())
        sites         list of sites: site 0 is the interstitial, then follows
                      the largest site, and then sites in decreasing order.
                      Each site is a list of tuples. Each tuple is the index
                      (i,j,k) into the map or grid.
        graph         NetworkX graph of the cells

        unit          physical units of various components
        P             (default) values of parameters

    Methods:

        map_sites(threshold)
                        label all sites, defined by the threshold. The threshold
                        value is stored with the object as the default. The default
                        can be explicitly set as P['threshold']
        save(filename)  save object.pickle
        load(filename)  restore object.pickle (or use d=Density(filename=<filename>))
        export()        write density to a file for visualization
        export_map()    write individual sites
    """

    # will probably break under multiple inheritance but I haven't figured out how to use super here
    _saved_attributes = Grid._saved_attributes + \
        ['map','sites','site_properties','equivalent_sites_index']

    # offsets to calculate first neighbours in the first octant
    __delta_first_octant__ = numpy.array([
        [0,0,1],
        [0,1,0],
        [1,0,0],
        [0,1,1],
        [1,0,1],
        [1,1,0],
        [1,1,1]])

    parameters_default = {'isDensity':False, 'threshold':None, 'MINsite':1}
    unit_default = {'length':'Angstrom', 'density':None}

    def __init__(self,grid=None,edges=None,filename=None,dxfile=None,
                 parameters=None,unit=None,metadata=None):
        """Adds information about sites to the grid. Sites are all
        cells with a density >= threshold.

        density = Density(kargs**)

        Sets up a Grid with additional data, namely the site map The
        threshold is given as key-value pair in the parameters
        dictionary and is assumed to be in the same units as the
        density.

        If the input grid is a histogram then it is transformed into a
        density.

        When neither grid nor edges are given then the density object
        can also be read from a pickled file (filename) or a OpenDX
        file (dxfile). In the latter case, care should be taken to
        properly set up the units and the isDensity parameter if the
        dx file is a density:

        >>> g = Density(dxfile='bulk.dx',parameters={'isDensity':True},
                    unit={'length':'Angstrom','density':'Angstrom'}, ....)
        """

        parameters = DefaultDict(self.parameters_default,parameters)
        if parameters['MINsite'] < 1:
            raise ValueError('MINsite must be > 0 (site must contain at least one cell)')
        if parameters['MINsite'] > 2:
            raise NotImplementedError("Currently only MINsite == 1 or 2 is supported.")
        unit = DefaultDict(self.unit_default,unit)
        metadata = DefaultDict({},metadata)

        # not necessary to initialize
        #self.graph = NX.Graph()      # nodes are indices (x,y,z) in isSite/grid
        #self.graph.name = 'density'  # graph.info() ...

        # All initialisations MUST come before calling Grid.__init__() because
        # of the crappy way that loading from a file is set up. (Otherwise I would overwrite
        # the loaded values.)
        self.map = None
        self.sites = None
        self.site_properties = None
        self.equivalent_sites_index = None

        Grid.__init__(self,grid=grid,edges=edges,filename=filename,dxfile=dxfile,
                      parameters=parameters,unit=unit,metadata=metadata)

        if not self.P['isDensity']:
            self.make_density()
        self.unit['threshold'] = self.unit['density'] # inconsistent with load() ?

    def map_sites(self,threshold=None):
        """Find regions of connected density and label them consecutively

        map_sites([threshold=<threshold>])


        threshold      Use the given threshold to generate the graph; the threshold
                       is assumed to be in the same units as the density.
                       (This updates the Density object's threshold value as well.)

        The interstitial has label '0', the largest connected subgraph
        has '1' etc. The sites (i.e.the list of indices into map/grid)
        can be accesed as Density.sites[label].
        """
        self.unit['threshold'] = self.unit['density']
        if threshold:
            self.P['threshold'] = threshold
        else:
            try:
                threshold =  self.P['threshold']
                if threshold is None:
                    raise ValueError
            except (KeyError,ValueError):
                raise ValueError("A threshold density is required to map the sites.")

        # first clean up
        try:
            self._remove_equivalence_sites()
        except ValueError:
            pass
        try:
            self._site_remove_bulk()
            warnings.warn("Removed bulk site in order to map sites. Add it again when this is done!")
        except (ValueError, TypeError):
            pass

        # map charts the sites. It starts out with the interstitial labeled
        # as 0 and everything else 1.
        self.map = numpy.where(self.grid >= self.P['threshold'],
                               1, SITELABEL['interstitial']).astype(numpy.int16)
        self._make_graph()
        self._label_connected_graphs()

    def map_hilo(self, lomin=0.0, lomax=0.5, himin=2.72):
        """**Experimental** mapping of low density sites together with high density ones.

        :Keywords:
          *lomin*
              low-density sites must have a density > *lomin* [0.0]
          *lomax*
              low-density sites must have a density < *lomax* [0.5]
          *himin*
              high-density sites must have a density > *himin* [2.72]
        """
        self.unit['threshold'] = self.unit['density']
        self.P['lomin'] = lomin
        self.P['lomax'] = lomax
        self.P['himin'] = himin
        self.P['threshold'] = self.P['himin']  # is this a good idea?

        # first clean up
        try:
            self._remove_equivalence_sites()
        except ValueError:
            pass
        try:
            self._site_remove_bulk()
            warnings.warn("Removed bulk site in order to map sites. Add it again when this is done!")
        except (ValueError, TypeError):
            pass

        # map charts the sites. It starts out with the interstitial labeled
        # as 0 and everything else 1.
        self.map = numpy.where(
            numpy.logical_or(
                numpy.logical_and(lomin < self.grid, self.grid < lomax), # low-density sites
                self.grid > himin),                                      # high-density sites
            1, SITELABEL['interstitial']).astype(numpy.int16)
        self._make_graph()
        self._label_connected_graphs()


    # 'site' functions (should make site a class)
    def site_occupancy(self,**labelargs):
        """Returns the labels and the average/stdev occupancy of the labeled site(s).

        labels, <N>, std(N) = site_coocupancy(include='all' | <int> | <list>)

        Average occupancy is the average number of water molecules on the site i:

           <N_i> = <n_i> * V_i

        where n_i is the average density of the site and V_i its volume.

        The label selection arguments are directly passed to
        site_labels() (see doc string).

        If the interstitial is included then 0,0 is returned for the
        interstitial site (so ignore those numbers).
        """
        labels = self.site_labels(**labelargs)
        props = self.site_properties.take(labels)
        return labels, props.occupancy_avg, props.occupancy_std

    def site_volume(self,**labelargs):
        """Returns the label(s) and volume(s) of the selected sites.

        labels, volumes = site_volume('all')

        The volume is calculated in the unit set in
        unit['length']. The label selection arguments are directly
        passed to site_labels() (see doc string).

        The volume of the interstitial (if included) is returned as 0
        (which is not correct but for technical reasons more
        convenient).
        """
        labels = self.site_labels(**labelargs)
        return labels, self.site_properties.take(labels).volume

    # methods to calculate site_properties
    def _site_occupancies(self,labels):
        factor = constants.get_conversion_factor('density',self.unit['density'],
                                                     self.unit['length'])
        Vcell = factor * numpy.linalg.det(self.delta)
        def occupancy(site):
            V = len(site) * Vcell
            if V == 0:  # interstitial
                return 0.0,0.0
            n = numpy.array([self.grid[s] for s in site])
            return numpy.average(n) * V, numpy.std(n) * V
        return numpy.array([occupancy(self.sites[l]) for l in labels])

    def _site_volumes(self,labels):
        Vcell = numpy.linalg.det(self.delta)
        return Vcell * numpy.array([len(self.sites[l]) for l in labels])

    def _site_centers(self,labels):
        """Returns geometric centers of sites or None for the interstitial."""
        m = self.midpoints
        def cellcenter(site):
            return numpy.array((m[0][site[0]],m[1][site[1]],m[2][site[2]]))
        def average_or_none(sites):
            try:
                return numpy.mean([cellcenter(site) for site in sites],axis=0)
            except ZeroDivisionError:
                return None
        return [average_or_none(self.sites[l]) for l in labels]

    def site_labels(self,include='default',exclude='default'):
        """Return a list of site labels, possibly filtered.

        L = site_labels(include=<inclusions>,exclude=<exclusions>)

        <inclusions> and <exclusions> consist of a list of site labels
        (integers) and/or keywords that describe a site selection. All entries
        in one list are logically ORed. All exclusions are then removed from
        the inclusions and the final list of site labels is returned as a numpy
        array. (As a special case, the argument need not be a list but can be a
        single keyword or site label).

        For convenience, some inclusions such as 'subsites' and
        'equivalencesites' automatically remove themselves from the exclusions.

        For standard use the defaults should do what you expect, i.e. only see
        the sites that are relevant or that have been mapped in a hopping
        trajectory.

        Set verbosity to 10 in order to see the parsed selection.

        <inclusions>
          'all'       all mapped sites, including bulk and subsites of
                      equivalent sites (but read the NOTE below: set exclude=None)
          'default'   all mapped sites, including bulk but excluding subsites
                      and interstitial
          'sites'     all mapped sites, excluding bulk and interstitial
                      (removes 'subsites' and 'equivalencesites' from
                      exclusions)
          'subsites'  all sites that have been renamed or aggreated into equivalence sites
          'equivalencesites'
                      only the equivalence sites
          int, list   site label(s)

        <exclusions>
          'default'    equivalent to ['interstitial','subsites']; always applied unless
                       exludsions=None is set!
          None         do not apply any exclusions
          'interstitial'
                       exclude interstitial (almost no reason to
                       ever include it)
          'subsites'   exclude sites that have been aggregated or
                       simply renamed as equivalence sites
          'equivalencesites'
                       exclude equivalence sites (and possibly include subsites)
          'bulk'       exclude the bulk site

        Provides the ordered list L of site labels, excluding sites listed in
        the exclude list. Site labels are integers, starting from '0' (the
        interstitial). These labels are the index into the site_properties[]
        and sites[] arrays.

        NOTE that by default the standard exclusions are already being applied
        to any 'include'; if one really wants all sites one has to set
        exclude=None.

        Exclusions are applied _after_ inclusions.

        'site' discards the bulk site, self.P['bulk_site']; this parameter is
        automatically set when adding the bulk site with site_insert_bulk().

        See find_equivalence_sites_with() for more on equivalence sites and
        subsites.
        """
        all_args = {'inclusions':
                        {None:['default'], 'default':['default'],
                         'all':['all'], 'sites':['sites'],
                         'subsites':['subsites'],
                         'equivalencesites':['equivalencesites']},
                    'exclusions':
                        {None:[], 'default':['interstitial','subsites'],
                         'interstitial':['interstitial'], 'bulk':['bulk'],
                         'subsites':['subsites'], 'equivalentsites':['equivalencesites'],
                         'equivalencesites':['equivalencesites']}
                    }
        def process(kind,arg):
            if type(arg) is not list:
                arg = [arg]
            _args = []
            for x in flatten(arg):
                if type(x) is int:          # add a site label
                    _args.append(x)
                elif x in all_args[kind]:   # add a keyword
                    _args.extend(all_args[kind][x])
                else:
                    raise ValueError(str(x)+" must be a (list of)  keyword(s) from "\
                                         +str(all_args[kind].keys())\
                                         +" or a (list of) integer(s) (i.e. site labels).")
            _args = numpy.unique(flatten(_args)).tolist()
            for i,x in enumerate(_args):
                try:   _args[i] = int(x)  # change unique-stringified site labels back to ints
                except ValueError:
                    pass
            return _args

        inclusions = process('inclusions',include)
        exclusions = process('exclusions',exclude)

        def remove_from(l,seq):                     # should be a set but that's not in py2.3 ... but should use 'from set import Set as set'
            """my_list = remove_from(my_list, [a,b,c,...])"""
            l = numpy.unique(l).tolist()            # must be unique and a list
            for x in seq:
                try:
                    l.remove(x)
                except ValueError:
                    pass
            return l

        ## TODO
        ## organize labels as set and make an sorted array at the end
        ##... site_labels = set(range(len(self.sites)))    # 'all'
        ##... site_labels.intersection_update(labels)
        ##... def add_labels(labels):
        ##        site_labels.union_update(labels)
        ##list is filtered below
        try:
            all_site_labels = range(len(self.sites))    # 'all', but list is filtered below
        except TypeError:
            return numpy.array([])                  # no sites to start with

        # TODO: horrible code: this should be cleaned up so that the
        # selections actually behave like proper boolean expressions
        # (intersection and union) and not an idiosyncratic mess
        site_labels = []
        for inclusion in inclusions:
            if inclusion == 'default' or inclusion is None:
                exclusions.append('subsites')
                site_labels.extend(all_site_labels)
            elif inclusion == 'all':
                site_labels.extend(all_site_labels)
            elif inclusion == 'sites':
                exclusions.extend(['bulk','interstitial'])
                exclusions = remove_from(exclusions,['subsites'])
                site_labels.extend(all_site_labels)
            elif inclusion == 'equivalencesites':
                exclusions = remove_from(exclusions,['equivalencesites'])
                site_labels.extend(self._labels_equivalencesites())
            elif inclusion == 'subsites':
                exclusions = remove_from(exclusions,['subsites'])
                site_labels.extend(self._labels_subsites())
            else:
                if inclusion not in all_site_labels:
                    raise ValueError('site label %d does not exist' % inclusion)
                site_labels.append(inclusion)  # site label

        inclusions = numpy.unique(inclusions)  # only for diagnostics
        exclusions = numpy.unique(exclusions)  # for the loop
        for exclusion in exclusions:
            if exclusion == 'interstitial':
                site_labels = remove_from(site_labels, [SITELABEL['interstitial']])
            elif exclusion == 'subsites':
                site_labels = remove_from(site_labels, self._labels_subsites())
            elif exclusion == 'equivalencesites':
                site_labels = remove_from(site_labels, self._labels_equivalencesites())
            elif exclusion == 'bulk':
                try:
                    site_labels = remove_from(site_labels, [self.P['bulk_site']])
                except KeyError:
                    pass
        msg(10,"site_labels(): using include = "+str(inclusions)+"  exclude = "+str(exclusions)+"\n")
        return numpy.unique(site_labels)

    def _labels_subsites(self):
        try:
            return self.site_properties.label[self.site_properties.equivalence_site != 0]
        except AttributeError:   # should not fail just because there are no equivsites
            return numpy.array([])

    def _labels_equivalencesites(self):
        try:
            return numpy.array(self.equivalent_sites_index.values())
        except AttributeError:  # should not fail just because there are no equivsites
            return numpy.array([])

    def subsites_of(self,equivsites,kind='sitelabel'):
        """Return subsites of given equivalence sites as a dict.

        dict <-- subsites_of(equivsites,kind='sitelabel')

        The dict is indexed by equivsite label. There is one list of subsites for
        each equivsitelabel.

        kind      'sitelabel':  equivsites are the sitelabels as uses internally; this is
                                the default because site_labels() returns these numbers and
                                so  one can  directly use the output from site_labels() as
                                input (see example)
                  'equivlabel': equivsites are treated as labels of equivalence sites;
                                these are integers N that typically start at 2
                  'name':       equivsites are treated as strings that are given as names
                                to sites; the default settings produce something like 'N*'

        EXAMPLES:

            dens.subsites_of(dens.site_labels('equivalencesites'))
            dens.subsites_of([2,5,10], kind='equivsites')
            dens.subsites_of('10*', kind='name')

        NOTE:
        * equivlabel == 0 is silently filtered (it is used as a merker for NO equivalence
          site)
        * empty equivalence sites show up as empty entries in the output dict; typically
          this means that one gave the wrong input or kind
        """
        _transform = {'sitelabel': self._sitelabel2equivlabel,
                      'equivlabel': asiterable,
                      'name': self._equivname2equivlabel,
                      }
        SP = self.site_properties
        subsites = {}
        try:
            equivsitelabels = _transform[kind](equivsites)
        except KeyError:
            raise ValueError('kind must be one of %r' % _transform.keys())
        # Must be a loop over equivlabels because one equivsite can contain many subsites
        # and I cannot do this as a direct index: the data structure is too crappy for that.
        for equivsite in equivsitelabels:
            if equivsite == 0:
                continue   # hack: we are using '0' as marker for NO equivsite so skip it
            subsites[equivsite] = SP[SP.equivalence_site == equivsite] # sites referencing the equivsite
        return subsites

    def _equivlabel2sitelabel(self,equivlabels):
        # XXX: not needed at the moment/untested
        return numpy.array([self.equivalent_sites_index[l] for l in asiterable(equivlabels)])

    def _equivlabel2equivname(self,equivlabels):
        """Returns the equivalence name strings for the equivalence labels."""
        return self.site_properties[[self.equivalent_sites_index[l] for l in
                                    asiterable(equivlabels)]].equivalence_name

    def _sitelabel2equivlabel(self,sitelabels):
        """Return equivalence labels corresponding to the internal sitelabels."""
        s = self.site_properties[asiterable(sitelabels)]
        return s[s.equivalence_label != 0].equivalence_label  # filter sites that are not equivsites

    def _equivname2sitelabel(self,equivnames):
        """Return internal sitelabel (169,170,..) corresponding to equivname ('2*', '3*')."""
        # XXX: not needed at the moment/untested
        SP = self.site_properties
        return numpy.ravel([SP[SP.equivalence_name == str(n)].label for n in asiterable(equivnames)])

    def _equivname2equivlabel(self,equivnames):
        """Return equivalence labels (2.3,..) corresponding to equivname ('2*', '3*')."""
        SP = self.site_properties
        return  numpy.ravel(
            [SP[SP.equivalence_name == str(n)].equivalence_label for n in asiterable(equivnames)])


    def site_insert_bulk(self,bulkdensity,bulklabel=SITELABEL['bulk'],force=False):
        """Insert a bulk site from a different density map as bulk site into this density.

        site_insert_bulk(bulkdensity)

        This is a bit of a hack. The idea is that one can use a site
        from a different map (computed from the same trajectory with
        the same grid!) and insert it into the current site map to
        define a different functional region. Typically, the bulk site
        is the largest site in bulkdensity (and has site label 1) but
        if this is not the case manually choose the appropriate
        bulklabel.

        The site is always inserted as the bulk site in the current density.

        Example:
        >>> bulkdensity = hop.interactive.make_density(psf,dcd,'bulk',delta=1.0,
                                atomselection='name OH2 and not within 4.0 of protein')
        >>> bulkdensity.map_sites(threshold=0.6)
        >>> density.site_insert_bulk(bulkdensity)
        >>> density.save()
        >>> del bulkdensity
        """
        # sanity checks
        if self.P.has_key('bulk_site') and not force:
            raise ValueError('The density already contains a bulk site. Use '
                               'force=True to override.')
        if self.map is None or bulkdensity.map is None:
            raise ValueError('Both densities must have had their map computed.')
        if self.map.shape != bulkdensity.map.shape:
            raise ValueError("The bulk density was defined on a different grid than this density.")
        if self.unit['threshold'] != bulkdensity.unit['threshold']:
            warnings.warn("The unit for the density (%s) is different from the unit "
                          "for the bulk density (%s).\n" %
                          (self.unit['threshold'], bulkdensity.unit['threshold']),
                          category=InconsistentDataWarning)
        # do the hack & update
        self.sites.insert(SITELABEL['bulk'],bulkdensity.sites[bulklabel])
        self._draw_map_from_sites()
        self._annotate_sites()
        self.P['bulk_site'] = SITELABEL['bulk']
        self.P['bulk_threshold'] = bulkdensity.P['threshold']

    def site_remove_bulk(self,force=False):
        """Cleanup bulk site."""
        self._site_remove_bulk(force=force)
        self._draw_map_from_sites()
        self._annotate_sites()

    def _site_remove_bulk(self,force=False):
        if not self.P.has_key('bulk_site') and not force:
            raise ValueError("No bulk site is recorded but force=True would do it")
        del self.sites[self.P['bulk_site']]
        del self.P['bulk_site']
        del self.P['bulk_threshold']

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

    def has_bulk(self):
        """Returns ``True`` if a bulk site has been inserted and ``False`` otherwise."""
        try:
            return self.P['bulk_site'] == SITELABEL['bulk']
        except KeyError:
            return False

    def masked_density(self,density,site_labels):
        """Returns only that portion of density that corresponds to
        sites; everything else is zeroed.

        masked = masked_density(density,sites)

        Arguments:

        density       a array commensurate with the map
        site_labels   label or list of site labels

        Results:

        Returns numpy array of same shape as input with non-site cells zeroed.
        """
        self._check_map()             # minimum sanity checks
        d = numpy.asarray(density)
        if not d.shape == self.map.shape:
            raise ValueError('Input "density" must have same shape as map.')
        try:
            site_labels[0]
        except (TypeError, IndexError): #  (numpy.int32 == scalar -> IndexError change in numpy...)
            site_labels = [site_labels]
        # logical OR all site maps and return corresponding entries in d
        return d * numpy.logical_or.reduce([self.map == label for label in site_labels])


    def export_map(self,labels='default',format='dx',directory=None,
                   value='density',combined=False,verbosity=3):
        """Write sites as a density file for visualization.

        export_map(**kwargs)

        labels='default'  Select the sites that should be exported. Can be
                           a list of numbers (site labels) or one of the keywords
                           recognized by site_labels() (qv). The interstitial is
                           always excluded.
        combined=False     True: write one file. False: write one file for each
                           site.
        format='dx'        Only dx format supported
        directory='site_maps'
                           Files are created in new directory, 'site_maps' by
                           default. File names are generated and indexed with
                           the label of the site. By default, 'site_maps' is
                           located in the same directory as the default filename.
        value= 'density'   Writes the actual density in the site.
               'threshold' The densities have the threshold value wherever the
                           site is defined. Note that the interstitial (label = 0)
                           is also written.
               <float>     Writes the value <float> into the site.
        verbosity=3        Set to 0 to disable status messages.

        Quick hack to write out sites. Each site can be written as a
        separate density file (combined=False) so that one can distinguish them
        easily in say VMD. Display with
           vmd site_maps/*.dx
        """
        set_verbosity(verbosity)
        self._check_map()        # minimum sanity checks
        threshold = self.P['threshold']
        # opt processing
        if directory is None:
            try:
                parentdir = os.path.dirname(self.filename())
            except ValueError:
                parentdir = os.path.curdir  # no default filename, so use cwd
            directory = os.path.join(parentdir,'site_maps')
        try:
            os.mkdir(directory)
        except os.error,e:
            if e.errno != errno.EEXIST:
                raise
        else:
            msg(3,'Created site map directory "%s".\n' % directory)
        # interstitial map MUST be excluded for 'painting over' to work
        site_labels = self.site_labels(labels,exclude=['interstitial'])
        if value == 'threshold':
            def sitemap(label,value=threshold):
                return numpy.where(self.map == label,value,0.0)
        elif value == 'density':
            def sitemap(label):
                return self.masked_density(self.grid,label)
        else:
            try:
                value = float(value)
            except ValueError:
                raise ValueError("value must be a number or 'threshold' or 'density'.")
            def sitemap(label,value=value):
                return numpy.where(self.map == label,value,0.0)

        def export_site_dx(g,filename,label,comments):
            components = dict(
                positions = OpenDX.gridpositions(1,g.shape,self.origin,self.delta),
                connections = OpenDX.gridconnections(2,g.shape),
                data = OpenDX.array(3,g),
                )
            dx = OpenDX.field('site map '+str(label),
                              components=components,comments=comments)
            dx.write(filename)
            msg(3,'Exported density %3s to file "%s".\n' % (str(label),filename))

        if combined:
            filename = os.path.join(directory,
                                    'sitemap_combined__%.2f.dx' % (threshold))
            g = sitemap(site_labels[0])      # first sitemap as basis (with its interstitial)
            for label in site_labels[1:]:
                smap = sitemap(label)        # set any values that are not interstitial in smap
                sites =  smap != SITELABEL['interstitial']
                g[sites] = smap[sites]       # 'paint over' previous values at sites
            comments = [
                'OpenDX density file written by export_map(combined=True),',
                '$Id$',
                '\t value at site is '+str(value),
                '\t threshold  = '+str(threshold),
                ]
            export_site_dx(g,filename,'all',comments)
        else:
            for label in site_labels:
                g = sitemap(label)
                filename = os.path.join(directory,
                                        'sitemap_%04d__%.2f.dx' % (label,threshold))
                comments = [
                    'OpenDX density file written by export_map(),',
                    '$Id$',
                    '\t value at site is '+str(value),
                    '\t threshold  = '+str(threshold),
                    '\t site_label = '+str(label),
                    ]
                export_site_dx(g,filename,label,comments)

    def find_equivalence_sites_with(self,reference,fmt='%d*',update_reference=True,
                                    use_ref_equivalencesites=False,
                                    verbosity=0,equivalence_graph='equivalence_graph.png'):
        """Find overlapping sites with a reference density and update site descriptions.

        Density.find_equivalence_sites_with(ref)

        :Arguments:
        ref         a Density object defined on the same grid
        fmt         python format string used for the equivalent_name, which should
                    contain %d for the reference label number (max 10 chars)
                    (but see below for magical use of xray water names)
        update_reference
                    True (default): Also update the site_properties in the
                    *reference* so that one can make graphs that highlight the
                    common sites. (This is recommended.)
                    False: don't change the reference
        use_ref_equivalencesites
                    True: use sites + equivalence sites from the reference density
                    False*: remove all equivalence sites als from the ref density
        verbosity   For verbosity >= 3 output some statistics; verbosity >=5 also
                    returns the equivalence graph for analysis; verbosity >= 7
                    displays the graph (and saves to equivalence_graph.png).

        An 'equivalence site' is a site that contains all sites that
        overlap in real space with another site in the reference
        density. This also means that two or more sites in one density
        can become considered equivalent if they both overlap with a
        larger site in the other density, and it is also possible that
        one creates 'equivalence' chains (0,a) <-> (1,b) <-> (0,c) <->
        (1,d) (although (0,a) ~<-> (1,d), and by construction (0,a) ~<->
        (0,c) and (1,b) ~<-> (1,d)), leading to extensive equivalence
        sites.

        When hopping properties are computed, an equivalence site is
        used instead of the individual sub sites.

        The equivalence sites themselves are constructed as new sites
        and added to the list of sites; their site numbers are
        constructed by adding to the total number of existing
        sites. Sub-sites are marked up by an entry of the equivalence
        site's site number in site_properties.equivalence_site.

        The common sites are consecutively numbered, starting at 2,
        from the one containing most sites to the one with fewest.

        The method updates Density.site_properties.equivalent_name
        with the new descriptor of the equivalent site. Equivalent
        site names are consecutively numbered, starting at 2, and can
        be optionally formatted with the fmt argument.

        However, if the reference density was built from an X-ray density AND if
        each site corresponds to single X-ray water molecule then the
        equivalence names contain the water identifiers eg 'W136' or
        'W20_W34_W36'.

        See the hop.sitemap.find_common_sites() function for more details.
        """
        original_verbosity = get_verbosity()  # change verbosity temporarily
        set_verbosity(verbosity)

        # Find sites in density 0 and 1 that overlap in space:
        # 1) get mapping
        # 2) interprete mapping as a graph with nodes the sites in each
        #    density; a node is a tuple (density,site), eg (0,23) or (1,117).
        # 3) mapping induces edges: the connected subgraphs constitute the common sites
        # 4) list connected subgraphs consecutively and use its ordinal as the common
        #    site label
        densities = [self, reference]                         # 0 = self, 1 = reference
        SELF,REF = 0,1

        def warn_and_remove(density):
            if density.equivalent_sites_index:
                warnings.warn('Density '+str(density)+' already contains equivalent sites, '
                              'which will be overwritten.',
                              category=OverwriteWarning)
                density.remove_equivalence_sites()
        warn_and_remove(self)
        if not use_ref_equivalencesites:
            warn_and_remove(reference)
        msg(5,"equivalence sites: use equivalence sites from reference: %r\n"
            % use_ref_equivalencesites)
        # if we use_ref_equivalencesites then we keep the reference map that has
        # equivalence sites painted in; thus all matching will be done against
        # equivalence sites which act as ordinary sites
        #
        # TODO: Still cannot keep equivalence_names for the NEW equiv sites.
        #       Perhaps a new naming scheme which concatenates site numbers or uses
        #       a concatentaion of all equiv names if they already exist?

        msg(5,"equivalence sites: finding mapping and analysing equivalence graph\n")
        m = find_common_sites(densities[SELF],densities[REF]) # mapping site_i(0) <--> site_k(1)
        edges = numpy.zeros((len(m),2,2),dtype=numpy.int16)   # (0,site_i(0)) <---> (1,site_k(1))
        edges[...,0] = [SELF,REF]       # identifier for density 0 and density 1
        edges[...,1] = m                # fill second field with corresponding node/site label
        ebunch = [map(tuple,e) for e in edges]  # make nodes hashable tuples
        g = NX.Graph()
        g.add_edges_from(ebunch)
        commonsites = [i for i in NX.connected_components(g)]  # each item: collection of equivalent sites

        # liz overlap
        overlap = find_overlap_coeff(densities[SELF],densities[REF])

        if update_reference:
            densities_to_update = {SELF:densities[SELF],REF:densities[REF]}
        else:
            densities_to_update = {SELF:densities[SELF]}

        # Book keeping: site_properties is the central and messy data
        # structure; see _annotate_sites()
        #
        # dirty update of sites and site_properties:
        #   1) add equivalence sites to self.sites
        #   2) self._draw_map_from_sites() --- overwrites sub-sites with equivalence
        #   3) self._annotate_sites() recreates self.site_properties with equivalence
        #      sites at end
        #   4) add additional markup for equivalence sites
        # NOTE: another _annotate_sites() will destroy the additional markup (perhaps
        #       put the additional markup into __annotate_sites(), too ?)
        # Note: The original sites do not know the actual label of the site (the index
        #       in site and site_propeties) but only the equivalence site label
        #       which, however, is identical across the densities. So add extra dict:
        for density in densities_to_update.values():
            density.equivalent_sites_index = dict()  # equiv. label --> index in sites

        # 1)
        msg(5,"equivalence sites: creating equivalence sites\n")
        for isite,commonsite in enumerate(commonsites):
            label = SITELABEL['bulk'] + isite + 1 # labels start after 'bulk'
            c = numpy.array(commonsite)           # array of sites [[0,23],[1,17],[0,2],...]
            for idensity,density in densities_to_update.items():
                sitelabels = c[c[:,0] == idensity][:,1] # [23,2,...] for idensity == SELF==0
                # would be simpler if sites was a numpy array...
                s = map(tuple,  numpy.concatenate(
                        [density.sites[l] for l in sitelabels.astype(int)]))
                density.sites.append(s)  # add compound site to list of sites
                sindex = len(density.sites) - 1 # index in sites and site_properties
                density.equivalent_sites_index[label] = sindex
        # 2) + 3) + 4)
        # (horrible hack to get xray waters into equivalence_name)
        msg(5,"equivalence sites: updating annotation\n")
        for idensity,density in densities_to_update.items():
            density._draw_map_from_sites() # density updated;equiv.sites overwrite subsites
            density._annotate_sites()      # recalculate all site stats
            for elabel,sindex in density.equivalent_sites_index.items():
                density.site_properties.equivalence_label[sindex] = elabel
                # construct a site name from x-ray waters in the reference density
                try:
                    isite = elabel - SITELABEL['bulk'] - 1   # XXX: argh, horrible...
                    c = numpy.array(commonsites[isite])      # subsites of elabel/isite
                    sref = c[c[:,0] == REF][:,1]             # only REF subsite labels
                    xraynames = reference.W(                 # only sites in REF...
                        reference.site2resid( sref ), format=True)  # .. add Wxxx identfiers
                    equivalence_name = "_".join(xraynames)
                except (AttributeError,NotImplementedError):
                    # default if we cannot find xray waters
                    equivalence_name = fmt % elabel
                    #raise
                density.site_properties.equivalence_name[sindex] = equivalence_name
            for isite,commonsite in enumerate(commonsites):
                label = SITELABEL['bulk'] + isite + 1 # labels start after 'bulk'
                c = numpy.array(commonsite)           # array of sites [[0,23],[1,17],[0,2],...]
                sitelabels = c[c[:,0] == idensity][:,1] # [23,2,...] for idensity == SELF==0
                density.site_properties.equivalence_site[sitelabels] = label

        # statistics and liz's stupid hack for probability overlap of equivalent sites
        if msg(3):
            msg(3,"equivalence sites: statistics for %d equivalent sites\n" \
                    % len(commonsites))
            msg(3,"self density sites are now labelled according to the remapped density: thus density.sites[n] will"
                  " NOT be equal to 'n' in 'sites' (0,n)")
            msg(3,"[%5s]  %5s  %5s %5s  |   %s   | %s | %s\n" % ('label','total','self','ref','sites','overlap coeff','total_overlap'))
            for isite,commonsite in enumerate(commonsites):
                label = SITELABEL['bulk'] + isite + 1    # labels start after 'bulk'
                labelstr = fmt % label
                c = numpy.array(commonsite)
                nsites = [0,0]
                for idensity in [SELF,REF]:
                    sitelabels = c[c[:,0] == idensity][:,1]
                    nsites[idensity] = len(sitelabels)
                # liz getting the probability overlap
                oc = overlap[isite]
                print oc
                msg(3,"[%5s]  %5d  %5d %5d  |   %s   | %s | %s\n" %
                    (labelstr,len(commonsite),nsites[SELF],nsites[REF],str(sorted(commonsite)),str(oc),str(sum(overlap))))
        if msg(7):
            msg(7,"Plotting the equivalent sites graph; blue in this density, red in reference density\n")
            import pylab
            n = numpy.array(g.nodes())
            node_color = numpy.where(n[:,0] == 0, 0.0,1.0)
            NX.draw_graphviz(g,node_size=200,node_color=node_color,alpha=0.2,font_size=7,linewidths=(0.1,),
                             labels=dict())
            NX.draw_graphviz(g,node_size=36,node_color=node_color,alpha=0.8,font_size=7)
            pylab.savefig(equivalence_graph)
        if msg(5):
            msg(5,"Returning equivalent-sites graph\n")
            set_verbosity(original_verbosity)
            return g
        set_verbosity(original_verbosity)

    def remove_equivalence_sites(self):
        """Delete equivalence sites and recompute site map."""
        self._remove_equivalence_sites()
        self._draw_map_from_sites()
        self._annotate_sites()

    def _remove_equivalence_sites(self):
        try:
            eqs_labels = self.site_labels(include='equivalencesites',exclude=None)
        except ValueError:
            return  # nothing to be done
        # Removing indexed items from a list is a bit iffy, and we can't do it sequentially.
        # This hack relies on the fact that the equivalence sites form one continuous block.
        try:
            first,last = eqs_labels[0], eqs_labels[-1]+1
            del self.sites[first:last]
            self.equivalent_sites_index = None
        except IndexError:
            return  # no eqs_labels

    def stats(self,data=None):
        """Statistics for the density (excludes bulk, interstitial, subsites).

        d = stats([data=dict])

        """
        if (not hasattr(self,'site_properties') or self.site_properties is None):
            raise AttributeError('Stats require site_properties annotation.')
        stats = dict()

        nodes = self.site_labels(exclude=['bulk','interstitial','subsites'])

        # General (meta data)
        try:
            stats['rho_cut'] = self.P['threshold']
            stats['rho_cut_bulk'] = self.P['bulk_threshold']
        except KeyError:
            warnings.warn("No bulk site defined", category=MissingDataWarning)
        stats['N_sites'] = len(nodes)
        stats['N_equivalence_sites'] = len(self.site_labels(include='equivalencesites',exclude=None))
        stats['N_subsites'] = len(self.site_labels(include='subsites',exclude=None))

        #------------------------------------------------------------
        # density stats
        #------------------------------------------------------------
        sp = self.site_properties[nodes]
        stats['site_volume_avg'] = numpy.average(sp.volume)
        stats['site_volume_std'] = numpy.std(sp.volume)
        stats['site_occupancy_rho_avg'] = numpy.average(sp.occupancy_avg)
        stats['site_occupancy_rho_std'] = numpy.std(sp.occupancy_avg)

        try:
            data['site_volume'] = sp.volume
            data['site_occupancy_rho_avg'] = sp.occupancy_avg
            data['site_occupancy_rho_std'] = sp.occupancy_std
        except TypeError:
            pass

        return stats

    def export3D(self,filename=None,site_labels='default'):
        """Export pdb and psf file of site centres for interactive visualization.

        >>> density.export3D()

        :Arguments:
        filename     prefix for output files:
                     <filename>.psf, <filename>.pdb, and <filename>.vmd
        site_labels  selects sites (See site_labels())

        The method writes a psf and a pdb file from the site map,
        suitable for visualization in, for instance, VMD. In addition,
        a VMD tcl file is produced. When it is sourced in VMD then the
        psf and pdb are loaded and site labels are shown next to the sites.

        Sites are represented as residues of resname 'NOD'; each site
        is marked by one 'ATOM' (of type CA) at the center of geometry
        of the site.

        Bulk and interstitial are always filtered from the list of
        sites because they do not have a well defined center.

        """
        if not hasattr(self,'site_properties') or self.site_properties is None:
            raise AttributeError('Requires site_properties ')

        site_labels = self.site_labels(site_labels,exclude=['interstitial','bulk'])

        self._write_psf(filename,site_labels)       # atoms are numbered consecutively...
        self._write_pdb(filename,site_labels)       # ..and residues correspond to sites
        self._write_vmd(filename,site_labels)       # tcl code for labels in VMD

    def _write_pdb(self,filename,site_labels):
        import Bio.PDB          # alternatively could use MDAnalysis, too, I guess
        B = Bio.PDB.StructureBuilder.StructureBuilder()
        B.init_structure('sites')
        B.init_model(0)
        B.init_seg('SITE')
        B.init_chain('A')
        # note that empty fields MUST be written as a blank ' ' (ie a single space)
        # or they are interpreted as altLoc specifiers named '' (weird...)
        # ATOM numbering is done consecutively here and in write_psf() and
        # there is only a single atom per residue.
        props = self.site_properties
        for node in site_labels:   # node is the label==resid and it must be an integer
            pos = props[node].center
            vol = props[node].volume
            occ = props[node].occupancy_avg
            commonlabel = props.equivalence_name[node].strip()
            if commonlabel:
                identity = 1.0
                aname = 'C'     # equiv/Common
            else:
                identity = 0.0
                aname = 'S'     # single site
            pdb_occupancy = occ    # this should be customizable and selected from
            pdb_beta = identity    # volume, occupancy, degree, identity
            B.init_residue('SIT',' ',node,' ') # choose same identifiers as in write_psf
            B.init_atom(aname,pos,pdb_beta,pdb_occupancy,' ','OH', element='O')
        io=Bio.PDB.PDBIO()
        s = B.get_structure()
        io.set_structure(s)
        pdbfile = self.filename(filename,'pdb')
        io.save(pdbfile)

    def _write_psf(self,filename,site_labels):
        """Pseudo psf with nodes as atoms and edges as bonds"""
        # Standard no CHEQ format for a Charmm PSF file:
        psf_ATOM_format = '%(iatom)8d %(segid)4s %(resid)-4d %(resname)4s '+\
                          '%(name)-4s %(type)4s %(charge)-14.6f%(mass)-14.4f%(imove)8d\n'

        psffilename = self.filename(filename,'psf')
        psf = open(psffilename,'w')
        psf.write('PSF\n\n')
        psf.write('%7d !NTITLE\n' % 2)
        psf.write('* Sites written by\n'+\
                  '* $Id$\n')
        psf.write('\n')

        # ATOMS
        psf.write('%6d !NATOM\n' % len(site_labels))
        segid = 'SITE'     # choose the same identifiers as in pdb
        resname = 'SIT'    # choose the same identifiers as in pdb
        atype = 'OH'       # choose the same identifiers as in pdb
        charge = 0
        mass = 1.0
        imove = 0            # no fixed 'atoms'
        props = self.site_properties
        for iatom,node in enumerate(site_labels):
            # atom numbering starts at 1, so iatom+1
            iatom += 1
            commonlabel = props.equivalence_name[node].strip()
            if commonlabel:
                aname = 'C'     # equiv/Common    same as in pdb
            else:
                aname = 'S'     # single site
            psf.write(psf_ATOM_format %
                      {'iatom':iatom, 'segid':segid, 'resid':node,
                       'resname':resname, 'name':aname, 'type':atype,
                       'charge':charge, 'mass':mass,'imove':imove} )
        # ignore all the other sections (don't make sense anyway)
        psf.close()

    def _write_vmd(self,filename,site_labels):
        vmdfilename = self.filename(filename,'vmd')
        PSF = self.filename(filename,'psf')    # regenerate filenames ....
        PDB = self.filename(filename,'pdb')    # ... quick and dirty
        vmd = open(vmdfilename,'w')
        vmd.write("mol new %(PSF)s type psf first 0 last -1 step 1 filebonds 1 autobonds 0 waitfor all\n"
                  "mol addfile %(PDB)s type pdb first 0 last -1 step 1 filebonds 0 autobonds 0 waitfor all\n"\
                  % vars())
        vmd.write("proc labelsites {} {\n")
        props = self.site_properties
        for node in site_labels:
            pos = props[node].center + numpy.array([0.5,0.5,0.5]) # shift a little bit
            label = props.equivalence_name[node].strip()    # if poossible, use a equivalence name
            if not label:
               label = str(props.label[node])               # .. otherwise use site label
            vmd.write('\tlappend sitelabels [draw text {%8.3f %8.3f %8.3f} "%s"]\n' % (pos[0],pos[1],pos[2], label))
        vmd.write("""}
proc delsitelabels {} {
\tforeach id $sitelabels {
\t\tdraw delete $id
\t}
}

mol delrep 0 top
mol representation CPK 1.000000 0.300000 8.000000 6.000000
mol color Name
mol selection {all}
mol material Opaque
mol addrep top

display projection orthographic
axes location off
color Display Background white

# draw labels
label textsize 1.0
draw color black
labelsites
puts "Labels can be deleted with 'delsitelabels'."
""")
        vmd.close()



    def _check_map(self):
        try:
            self.map[0]
            if not self.map.shape == self.grid.shape:
                raise ValueError
        except (TypeError,AttributeError,ValueError):
            raise ValueError('The map has not been generated. Run map_sites() first.')

    def _make_graph(self):
        """Connect nodes that are in each others' first shell. Excludes cells that
        are not connected to any other cells, i.e. isolated cells will be counted
        as interstitial if MINsite > 1.

        make_graph()

        Populates Density.graph with nodes (cell coordinates) and
        edges (tuples of nodes).
        """
        # Currently excludes cells that are not connected to any other
        # cells. THIS MAY HAVE TO BE CHANGED (add [0,0,0] to __delta_ and
        # also check export_map()).

        self.graph = NX.Graph()      # nodes are indices (x,y,z) in map/grid
        self.graph.name = 'density'  # graph.info() ...

        # site index list:
        # [[x1,x2,...],[y1,y2,...],[z1,z2,...]] array of indices -->
        # belonging to sites (x1,y1,z1), (...)
        sidx = map(tuple,
                   numpy.array(numpy.where(self.grid >= self.P['threshold'])).T)
        sidx.sort()          # comes sorted but just make sure (important!)
        for site in sidx:    # TODO: Optimize these loops!
            for neighbour in self._shell(site):
                try:
                    if self.map[neighbour]:
                        self.graph.add_edge(site,neighbour)
                except IndexError:
                    pass     # neighbour outside grid (probably should pad map)
        if self.P['MINsite'] == 1:
            self.graph.add_nodes_from(sidx)  # include ALL high-density regions


    def _shell(self,site):
        """list of indices of neighbour cells in the (+,+,+) octant

        This will exhaust all cells because the site index list is
        sorted in the same direction, ie we are looking for edges in
        the 'forward' direction.

        Note:
        * The site itself is excluded from the neighbours.
        * Periodic boundary conditions have NOT been taken care off
        """
        return map(tuple, site + self.__delta_first_octant__)  # site must be a tuple

    def _label_connected_graphs(self):

        """finds all connected subgraphs, sorts the m by size (largest
        first), and labels the cells in the map by the rank of the
        corresponding subgraph

        The interstitial has label '0', the largest connected subgraph
        has '1' etc. The sites (i.e.the list of indices into map/grid)
        can be accesed as Density.site[label].

        Note that ISOLATED GRID CELLS do not show up as connected
        components and are also zeroed out in the map: they are
        COUNTED AS INTERSTITIAL. (One could change this by initially
        marking up the map with -1 and then looking for the -1 at the
        end of this function.)
        """

        self.sites = []
        components = list(NX.connected_components(self.graph))  # this does the hard work
        for component in components:
            self.sites.append(list(component))
        volume_list=map(len,self.sites) #pull all of the 'volumes'
        volume_sites_list=list(izip(volume_list,self.sites)) #crudely zip the volumes with the site lists for sorting
        volume_sites_list.sort() #should sort by the first element, the "volume"
        volume_sites_list.reverse()
        self.sites = [site[1] for site in volume_sites_list] #hopefully returns the sorted sites
        self.sites.insert(SITELABEL['interstitial'],[])   # placeholder for interstitial
        self._draw_map_from_sites()
        self._annotate_sites()


    def _draw_map_from_sites(self):
        """Label cells in the map based on the site list.

        _draw_map_from_sites()

        Note that later sites (higher label) overwrite earlier
        sites. This is only important if the site list was manually
        manipulated and some sites actually overlap.
        """
        # start with clean map (if MINsite>1 then isolated density becomes interstitial
        self.map = SITELABEL['interstitial'] * numpy.ones(self.grid.shape,dtype=numpy.int16)
        for label,sidx in enumerate(self.sites):
            for s in sidx:
                self.map[s] = label

    def _annotate_sites(self):
        """List of properties associated with each site. Overwrites previous site_properties."""
        # make sure that this also includes equivalence sites and their sub-sites
        labels = self.site_labels(include='all',exclude=None)
        occ_avg,occ_std = self._site_occupancies(labels).T[[0,1]]
        equivalence_label = [0] * len(labels)
        equivalence_site = [0] * len(labels)       # link to 'equivalence_site', see find_equivalence_sites_with()
        ##equivalence_name = [' '*10] * len(labels)  # fieldwidth 10
        equivalence_name = [None] * len(labels)    # want any-length strings
        empty_centers = [None] * len(labels)       # workaround to get an 'object' record
        # Define and initialize the important data structure 'site_properties'.
        # See also find_equivalent_sites_with()
        self.site_properties = numpy.rec.fromarrays([
            labels,                      # site number (id)
            self._site_volumes(labels),  # volume
            occ_avg,occ_std,             # occupancy: avg,stdev
            empty_centers,               # geometric centre
            equivalence_label,           # new label given to an equivalence site
            equivalence_site,            # label (int) of the equivalence site that contains this site
            equivalence_name,            # string representation of the label
            ],
            names='label,volume,occupancy_avg,occupancy_std,center,equivalence_label,equivalence_site,equivalence_name')
        self.site_properties.center[:] = self._site_centers(labels)  # now fill with 3-arrays
        self.site_properties.equivalence_name[:] = ''

    def __repr__(self):
        features = [str(self.grid.shape)+' bins',]  # always available
        try:
            fn = 'default filename "'+str(self.filename())+'"'
            features.append(fn)
        except ValueError:
            pass
        try:
            nsites = len(self.sites) - 1
            features.append(str(nsites)+' sites')
        except TypeError:
            pass
        try:
            nequivsites = len(self.equivalent_sites_index)
            features.append('including '+str(nequivsites)+' equivalence sites')
        except TypeError:
            pass
        return '<hop.sitemap.Density density with '+', '.join(features)+'>'


def remap_density(density,ref,verbosity=0):
    """Transform a Density object to a grid given by a reference Density.

    >>> newdensity = remap_density(old,ref)

    The user is repsonsible to guarantee that:
    * the grid spacing is the same in both densities
    * the grids only differ by a translation, not a rotation

    :Arguments:
    old          Density object with site map
    ref          reference Density object that provides the new grid shape
    verbosity=0  increase to up to 3 for status and diagnostic messages

    :Returns:
    newdensity   Density object with old's density and site map transformed
                 to ref's coordinate system. It is now possible to manipulate
                 newdensity's and ref's arrays (grid and map) together, e.g.

                 >>> common = (newdensity.map > 1 & ref.map > 1)
                 >>> pairs = newdensity.map[common], ref.map[common]

    Note that this function is not well implemented at the moment and
    can take a considerable amount of time on bigger grids
    (100x100x100 take about 3 Min).

    An implicit assumption is that the two coordinate systems for the
    two grids are parallel and are only offset by a translation. This
    cannot be checked based on the available data and must be
    guaranteed by the user. RMS-fitting the trajectories is sufficient
    for this to hold.

    BUGS:

    * This is not a good way to do the remapping: It requires parallel
      coordinate systems and the exact same delta.
    * It is slow.
    * It would be much better to interpolate density on the reference grid,
    """
    set_verbosity(verbosity)  # set to 0 for no messages

    try:
        ref.map
        ref.grid
        density.map
        density.grid
    except AttributeError:
        raise TypeError('Both density and ref need to be proper Density objects')

    if not numpy.all(numpy.abs(ref.delta - density.delta) < 1e-4):
        warnings.warn('The grid spacings are not the same; this is probably WRONG.',
                      category=InconsistentDataWarning)

    D=numpy.rank(ref.map)
    newgrid = numpy.zeros(ref.grid.shape)
    newmap = SITELABEL['interstitial'] * numpy.ones(ref.map.shape, dtype=numpy.int16)
    # lookup table for index transformation: i,j,k |--> i',j',k' = t(i,j,k)
    # (instead of a linear transformation; also makes it easier to flag outliers)
    t_table = [numpy.digitize(density.midpoints[axis],ref.edges[axis])-1 \
                            for axis in xrange(D)]
    for axis in xrange(D):    # mark outliers with -1 (and filter in transformed())
        x = t_table[axis]
        x[ (x<0)|(x>=ref.grid.shape[axis]) ] = -1

    def transformed(ijk=None):
        """Returns transformed index triplett or None if out of bounds or no input.
        >>> i',j',k' = transformed((3,4,5))
        """
        if ijk is None:
            return None     # called on empty site
        try:
            u,v,w = t_table[0][ijk[0]],t_table[1][ijk[1]],t_table[2][ijk[2]]
        except IndexError:
            return None     # outside the table
        if u == -1 or v == -1 or w == -1:
            # only return indices that are in bounds (outliers == -1)
            return None
        return u,v,w

    # Remap the density (moderately fast but not great):
    # pre-compute pairs (ijk,t(ijk)), only if t(ijk) != None then assign all at once
    msg(3,"Transforming %d indices...\n" % numpy.product(ref.grid.shape))
    # TODO OPTIMIZE:
    # next line is a bottleneck and takes long (the if-part is not the problem)
    idx = numpy.array([(ijk, transformed(ijk)) \
            for ijk in numpy.ndindex(density.map.shape) if transformed(ijk)])
    msg(3,"Remapping the density old --> new...\n")
    newgrid[idx[:,1,0],idx[:,1,1],idx[:,1,2]] = density.grid[idx[:,0,0],idx[:,0,1],idx[:,0,2]]

    msg(3,"Remapping the site list...\n")
    newsites = [unique_tuplelist([transformed(ijk) for ijk in site]) for site in density.sites]
    #DensityClass = type(density)     # can be Density or PDBDensity
    DensityClass = hop.sitemap.__getattribute__(type(density).__name__)
    newdensity = DensityClass(grid=newgrid,edges=ref.edges[:],
                         parameters=density.P.copy(),
                         unit=density.unit.copy(),metadata=density.metadata.copy())
    for attr in ('_xray2psf','_psf2xray'):
        try:
            # add xray2psf translation table if it is a PDBDensity
            newdensity.__setattr__(attr,density.__getattribute__(attr).copy())
        except AttributeError:
            pass
    newdensity.sites = newsites       # hack in new sites instead of calc. from density
    newdensity._draw_map_from_sites() # build map from transformed sites
    newdensity._annotate_sites()      # recalculate site props, esp. centers
    newdensity.P['remapped_source'] = str(density)
    newdensity.P['remapped_reference'] = str(ref)
    return newdensity

def unique_tuplelist(x):
    """Sort a list of tuples and remove all values None"""
    if len(x) == 0:
        return []
    tmp = numpy.empty(len(x),dtype=numpy.object_)  # must pre-allocate to keep dtype
    tmp[:] = x                                     # array with tuples as entries
    tmp.sort()
    if tmp[0] is None:          # None sorts into first position
        takefirst = False       # hack to eliminate None from result
    else:
        takefirst = True
    idx = numpy.concatenate(([takefirst],tmp[1:]!=tmp[:-1])) # remove duplicates
    return tmp[idx].tolist()    # here we want a list

def find_common_sites(a,b,use_equivalencesites=None):
    """Find sites that overlap in space in Density a and b.

    m = find_common_sites(a,b)

    :Arguments:
    a        Density instance
    b        Density instance

    :Returns:
    array of mappings between sites in a and b that overlap
    m[:,0]     site labels in a
    m[:,1]     site labels in b
    dict(m)    translates labels in a to labels in b
    dict(m[:,[1,0]])
               translates labels in b to labels in a

    """
    try:
        if a.map.shape != b.map.shape:
            raise ValueError('a and b are not defined on the same grid: use remap_density()')
        if not numpy.all([numpy.all(a.edges[axis] == b.edges[axis]) \
                              for axis in xrange(numpy.rank(a.map))]):
            raise ValueError('a and b do not superimpose in space (different edges)')
    except AttributeError:
        raise TypeError('Both densities need to be proper Density objects with map and edges attributes.')

    # make a list of sites that appear in the same position and are not bulk
    common = (a.map > SITELABEL['bulk']) & (b.map > SITELABEL['bulk'])
    m = numpy.array([a.map[common],b.map[common]]) # multiple entries == overlap vol
    m = numpy.unique(map(tuple,m.transpose()))     # remove multiple entries
    if not numpy.any(m):
        m = numpy.array([[],[]])                   # return empty mapping
    return m

def find_overlap_coeff(a,b):
    """Find sites that overlap in space in Density a and b.

    m = find_overlap_coeff(a,b)

    :Arguments:
    a        Density instance
    b        Density instance

    :Returns:
    array sites in a and b that overlap
       and array of probability of overlap for overlapped sites
    m[:,0]     site labels in a
    m[:,1]     site labels in b
    oc         amount of overlap
    """
    #### liz hack to beautiful code
    common = (a.map > SITELABEL['bulk']) & (b.map > SITELABEL['bulk'])
    m = numpy.array([a.map[common],b.map[common]]) # multiple entries == overlap vol
    m = numpy.unique(map(tuple,m.transpose()))     # remove multiple entries
    oc = numpy.zeros(len(m[:,0]))
    for isite, i in enumerate(m):
        coeff = 0
        #a_cellvol = a.delta[0][0]*a.delta[1][1]*a.delta[2][2]  # calc vol of grid cell for a
        #b_cellvol = b.delta[0][0]*b.delta[1][1]*b.delta[2][2]  # calc vol of grid cell for b
        sum_a = a.grid[a.map > SITELABEL['bulk']].sum()
        sum_b = b.grid[b.map > SITELABEL['bulk']].sum()
        #f_a = constants.get_conversion_factor('density', a.unit['length'], a.unit['density'])
        #f_b = constants.get_conversion_factor('density', b.unit['length'], b.unit['density'])
        #a_bulk_density = (a.site_occupancy()[1][0]/a.site_volume()[1][0]) * f_a
        #b_bulk_density = (b.site_occupancy()[1][0]/b.site_volume()[1][0]) * f_b
        #print a_cellvol, b_cellvol
        for j in a.sites[i[0]]:
                for k in b.sites[i[1]]:
                        if j == k:
                                #print a.grid[j], b.grid[k],j,i[0],i[1], a_bulk_density, b_bulk_density
                                if a.grid[j] < b.grid[k]:
                                        # density of grid cell point normalized to bulk density in a single grid cell
                                        coeff = coeff + (a.grid[j]/sum_a) #/a_bulk_density)
                                else:
                                        coeff = coeff + (b.grid[k]/sum_b) #b_bulk_density)
        oc[isite] = coeff
    return oc

