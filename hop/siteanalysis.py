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

"""
Advanced analysis of site properties --- :mod:`hop.siteanalysis`
================================================================

The module provides a framework on which to build analysis tasks on a
per-site basis.

.. warning:: Experimental code.

"""
from __future__ import absolute_imports

import warnings

import numpy

from . import MissingDataError
from . import trajectory
from . import utilities
from .utilities import set_verbosity, msg, iterable, IntrospectiveDict
from .constants import SITELABEL


## For launching into the debugger, add 'debug_here()' in the code
#from IPython.Debugger import Tracer; debug_here = Tracer()

class SiteArray(utilities.Saveable):
    """Base class to implement containers that hold observables from a trajectory.

    The base class is not useful on its own and must be derived from.

    :API:

    Subclasses should define at least the following so that the plot functions work:

      self._data         #  numpy.array with sites on axis=0 (including 'rubbish' rows etc)
      self.data          #  typically a view on _data; the data that makes sense to the user
      self.distribution  # a 'normalized' version of the data (probably only sensible for histograms)

    See comments in __init__ about self.site2indx on indexing in self._data.

    Subclasses must call the super init at the start of their own __init__():

        super(cls,self).__init__(name,sites,maxsites,parameters,properties=properties,**kwargs)

    The following methods must be defined:

       add()        Add data from a time step.
       reset()      Clear all data and be ready to acquire new data with the next add() call.
    """

    _saved_attributes = 'all'
    default_plotproperties = {'xlabel':"x", 'ylabel':None,
                              'title': "Site-resolved distributions",
                              'legend_fontsize':6,
                              'xticks_indexed':False,
                              'im_yticks_site':True, 'im_yticks_fontsize':6 # imshow
                              }  # can be changed by subclasses and properties argument
    fixed_plotproperties = {}    # cannot be changed by subclasses

    def __init__(self,*args,**kwargs):
        ## super(SiteArray,self).__init__(*args,**kwargs) # TODO: Saveable does not work this way yet!
        self.name = args[0]          # TODO: change this crap to keywords
        self.sites = args[1]
        self.maxsites = args[2]
        self.parameters = args[3]

        self.Nsites = Nsites = len(self.sites)

        # Setting up the translation table 'site --> index in array':
        #   self.site2indx[]
        # Lookup table for the index in the results data array corresponding to a
        # site label (=the number stored in the hopping trajectory and the hopping graph)
        # IMPORTANT: site labels that are NOT in <sites> point to the index 'Nsites' and
        # we add one additional data row = Nsites into which all data are dumped
        # that are not selected; this allows for more efficient, vectorized array access;
        # see longer comment in SiteHistogramArray.__init__().
        site2indx = Nsites * numpy.ones(self.maxsites+1, dtype=int)  # initialize to Nsites
        site2indx[self.sites] = numpy.arange(Nsites)            # set selected sites
        self.site2indx = site2indx    # maps site label to its array row (=indx)
        self.indx2site = self.sites   # <-- without uncollected site array aindx=Nsites

        # properties: customization of plots etc
        self.properties = self.default_plotproperties.copy()
        arg_properties = kwargs.pop('properties',{})
        self.properties.update(arg_properties)
        self.properties.update(self.fixed_plotproperties)

        # specialized setup must be done in subclasses, in particular self._data
        pass

    def mean(self):
        """Mean for each site."""
        return self.data.mean(axis=1)

    def std(self):
        """Standard deviation for each site."""
        return self.data.std(axis=1)

    def mean_std(self):
        "Returns site labels, mean, and standard deviation for each site."""
        return self.sites, self.mean(), self.std()

    def plot(self,filename=None,sites=None,normalize=True,format='pdf',with_legend=True,
             use_midpoints=True,step=None,**kwargs):
        """Plot the distribution(s) or histogram(s) for sites.

        :Arguments:
        filename      write to file filename; suffix determines type
        sites         site label or list of site labels; if None then all
                      collected sites are plotted.
        normalize     True: distribution; False: histogram
        format        determine output format if filename is derived from the
                      default filename
        with_legend   True | False
        use_midpoints True: plot distribution at midpoints of histogram bins
                      False: plot at lower edge
        step          plot only every <step> datapoints (along x); None: auto

        :Limitations:
        * sites must be a list of integers and a subset of Collection.sites.
          It would be nice if one could use the same selection syntax as for
          the constructor (which is just passed on to
          Density.site_labels()). The reason why this is not implemented in
          plot is because when loading from a pickled file the density is not
          available and neither is the site_labels() method.
        * The legend overflows with many sites; it should also show the
          equivalence site name if there is one.
        """
        import pylab
        if normalize:
            data = self.distribution
        else:
            data = self.data
        # sanity check
        try:
            data[0]
            if len(self.site2indx) == 0: raise
            if len(self.sites) == 0: raise
        except:
            if self.distribution is None and self.data is not None:
                raise ValueError("No normalized data available, try normalize=False.")
            else:
                raise MissingDataError("Cannot plot histograms because data are missing. Run compute() first.")

        step = self.autoset_step(step)
        data = data[:,::step]

        if sites is None:
            sites = self.sites
        elif not iterable(sites):
            sites = [sites]
        pylab.clf()
        if use_midpoints:
            X = self.midpoints
        else:
            X = self.edges[:-1]
        try:
            lines = [pylab.plot(X, data[self.site2indx[s]], label=str(s))
                     for s in sites]
        except IndexError:
            missing_sites = [s for s in sites if s not in self.sites]
            raise ValueError("No data for sites "+str(missing_sites)+
                             " (probably they need to be included in SiteAnalysis(sites=...) first.")
        pylab.xlabel(self.properties['xlabel'])
        pylab.title(self.properties['title'])
        if with_legend:
            pylab.legend(loc='best', prop=pylab.matplotlib.font_manager.FontProperties(
                    size=self.properties['legend_fontsize']))
        try:
            fig = self.filename(filename,ext=format,use_my_ext=True)
            pylab.savefig(fig)
            msg(3,"Saved figure to %s.\n" % fig)
        except ValueError:
            msg(1,"Cannot save figure without a filename.\n")

    def imshow(self,filename=None,sites=None,normalize=True,format='pdf',step=None,**kwargs):
        """Plot the distribution(s) or histogram(s) for sites as color grid.

        :Arguments:
        filename      write to file filename; suffix determines type
        sites         site label or list of site labels; if None then all
                      collected sites are plotted.
        normalize     True: distribution; False: histogram
        format        determine output format if filename is derived from the
                      default filename
        step          plot only every <step> datapoints (along x); None: auto
        **kwargs      keyword arguments are passed directly to pylab.imshow()

        :Limitations:
        * sites must be a list of integers and a subset of Collection.sites.
          It would be nice if one could use the same selection syntax as for
          the constructor (which is just passe on to
          Density.site_labels()). The reason why this is not implemented in
          plot because when loading from a pickled file the density is not
          available and neither is the site_labels() method.
        """
        import pylab
        # kwargs are passed to pylab.imshow():
        kwargs.setdefault('aspect','auto')
        kwargs.setdefault('interpolation','nearest')
        if normalize:
            data = self.distribution
            kwargs.setdefault('vmin',0.0)
            kwargs.setdefault('vmax',1.0)
        else:
            data = self.data
        # sanity check
        try:
            data[0]
            if len(self.site2indx) == 0: raise
            if len(self.sites) == 0: raise
        except:
            raise MissingDataError("Cannot plot histograms because data are missing. Run compute() first.")
        step = self.autoset_step(step)
        data = data[:,::step]

        # sites list: note that sites[row] is the proper label of this row in the image
        if sites is None:
            sites = self.sites
        elif not iterable(sites):
            sites = numpy.array([sites])
        else:
            sites = numpy.asarray(sites)
        xdelta = self.edges[1] - self.edges[0]
        shifted_edges = self.edges - 0.5*xdelta
        pylab.clf()
        try:
            img = pylab.imshow(data[[self.site2indx[s] for s in sites]],
                               extent=(shifted_edges[0],shifted_edges[-1],len(sites),0),
                               origin='upper',
                               **kwargs)
        except IndexError:
            missing_sites = [s for s in sites if s not in self.sites]
            raise ValueError("No data for sites "+str(missing_sites)+
                             " (probably they need to be included in SiteAnalysis(sites=...) first.")
        ax = img.get_axes()
        if self.properties['xticks_indexed']:
            xintegerLocator = pylab.matplotlib.ticker.IndexLocator(base=1,offset=-0.5*xdelta)
            ax.xaxis.set_major_locator(xintegerLocator)
        if self.properties['im_yticks_site']:
            _sites_as_str=sites.astype('|S10')
            def _idx2site(x,pos):
                try:
                    return _sites_as_str[int(x)]
                except IndexError:
                    return ''
            idx2siteFormatter = pylab.matplotlib.ticker.FuncFormatter(_idx2site)
            yintegerLocator = pylab.matplotlib.ticker.IndexLocator(base=1,offset=-0.5)
            ax.yaxis.set_major_locator(yintegerLocator)  # have sep. locs for diff axes (y and x)!
            ax.yaxis.set_major_formatter(idx2siteFormatter)
            for label in ax.yaxis.get_ticklabels():
                label.set_fontsize(self.properties['im_yticks_fontsize'])
            pylab.ylabel('site label')
        else:
            pylab.ylabel(self.properties['ylabel'])
        pylab.colorbar()
        pylab.xlabel(self.properties['xlabel'])
        pylab.title(self.properties['title'])
        pylab.draw_if_interactive()
        try:
            fig = self.filename(filename,ext=format,use_my_ext=True)
            pylab.savefig(fig)
            msg(3,"Saved figure to %s.\n" % fig)
        except ValueError:
            msg(1,"Cannot save figure without a filename.\n")

    def autoset_step(self,step):
        if step is not None:
            return step
        maxwidth = 800               # plot graphs up to <maxwidth> pixels wide
        ndata = len(self.data)
        if ndata > maxwidth:         # autoset step
            step = ndata/maxwidth
            msg(4,'Auto-set step to %(step)d.\n' % locals())
        else:
            step = 1
        return step


class SiteTrajectoryArray(SiteArray):
    """The whole trajectory, split by site.

    Organized as one BIG numpy array. The first axis correponds to the site
    index, the second axis corresponds to the frames. Note that the array is
    kept in memory. If the trajectory is long and/or many atoms have been
    selected then the machine may crash due to insufficient memory.
    """

    fixed_plotproperties = {'xlabel':r"time $t$/ps",
                            'xticks_indexed':False,
                            }
    _excluded_attributes = ['_next_frame_number']   # hack for Saveable

    def __init__(self,name,sites,maxsites,parameters,properties=None,**kwargs):
        """Initialize the SiteTrajectoryArray.

        A = SiteTrajectoryArray(sites,maxsites,properties=dict(numframes=<int>))

        :Arguments:
        name           name of the observable
        sites          list of all site labels
        maxsites       highest site label encountered in the trajectory
        parameters     dict(numframes=<number of frames that are going to be saved>,
                            totaltime=<length of trajectory>,
                            delta_t=<length of one frame>)
        properties     dict to customize the plot (xlabel, legend_fontsize, ...).
        """
        super(SiteTrajectoryArray,self).__init__(name,sites,maxsites,parameters,
                                                 properties=properties,
                                                 **kwargs)
        self.numframes = self.parameters['numframes']
        self.totaltime = self.parameters['totaltime']
        self.delta     = self.parameters['delta_t']
        self.edges = self.delta * numpy.arange(self.numframes+1)   # only edges makes really sense
        self.midpoints = 0.5*(self.edges[1:] + self.edges[:-1])    # but for compatibility...

        # allocate array (with one additional row to dump data from not-selected sites)
        msg(4,'SiteTrajectoryArray.__init__(): '
            'Allocating float data array of size %d x %d [%.1f MiB]' % \
                (self.Nsites+1,self.numframes,(self.Nsites+1)*self.numframes*32.0/(1024**2)))
        self._data = numpy.zeros((self.Nsites+1,self.numframes),dtype=float)

        # keep track of frame number internally
        # TODO: This makes save() fail (need to filter type('method-wrapper'))
        # (put this into _excluded_attributes[])
        self._next_frame_number = self._frame_number_generator().next

    def _frame_number_generator(self):
        """Supply a new frame number for each call of add()."""
        for i in xrange(self.numframes):
            yield i

    def add(self,site_data,sites):
        """Add one frame of data to the site trajectory.

        site_data     numpy array; axis=0 indexes the sites
                      For scalars, site_data.shape == (N,), with N == len(sites)
                      For vectors, site_data.shape == (N,3).
        sites         numpy array of the site labels that correspond to the slices
                      along axis=0 in site_data
        """
        assert len(sites) == len(site_data), "incompatible sites and site_data"
        try:
            # self.site2indx[sites] adds any data for not-selected sites to the 'rubbish' row
            self._data[self.site2indx[sites],self._next_frame_number()] = site_data
        except StopIteration:
            raise RuntimeError("Trying to read more frames than specified in numframes=%d. "
                               "You need to 'reset()' the observable '%s'."
                               % (self.numframes,self.name))
        #debug_here()

    def reset(self):
        """Erases all accumulated data and restarts the frame counter."""
        self._data *= 0
        self._next_frame_number = self._frame_number_generator().next

    def autocorrelation(self,start=None,stop=None,step=None,**kwargs):
        """Calculates the auto correlation function for all site trajectories."""
        from .utilities import autocorrelation_fft as ACF
        return numpy.array([ACF(traj,**kwargs)
                              for traj in self.data[:,start:stop:step]])

    def averaged_autocorrelation(self,step=None,**kwargs):
        """Calculates the ACF or each site by resampling from the whole trajectory.

        mean(acf), standardev(acf) = averaged_autocorrelation(**kwargs)

        :Arguments:
        step            only take every <step> from the trajectory (None == 1)
                        ??? step > 1 seems to take LONGER ???
        length          length (in frames) of the ACF (default: 1/2*len(series))
        sliding_window  repeat ACF calculation every N frames (default: len(series)/100)

        :Returns:
        mean_acf        average over all resampled acfs per site,  shape = (Nsites,length)
        std_acf         standard deviation or the resampled acfs,  shape = (Nsites,length)

        See also for kwargs:
        """
        from .utilities import averaged_autocorrelation as avACF
        tmp = numpy.array([avACF(traj,**kwargs) for traj in self.data[:,::step]])  # [(mean,std),...]
        return tmp[:,0], tmp[:,1]

    #averaged_autocorrelation.func_doc += utilities.averaged_autocorrelation.func_doc

    def _calculate_normalization(self):
        """Normalize each site to its highest value."""
        self.normalization = numpy.max(self.trajectory,axis=1) # does not include rubbish rows
        self.normalization[self.normalization == 0] = -1  # fix for normalization of empty rows

    def distribution():
        doc = """Numpy array of the trajectory, with each site normalized by its maximum value."""
        # implement as computed attribute to save memory for large trajectories
        def fget(self):
            try:
                len(self.normalization)
            except:
                self._calculate_normalization()
            return self.trajectory/self.normalization[...,numpy.newaxis]
        return locals()
    distribution = property(**distribution())

    def trajectory():
        doc = """Numpy array with N sites on axis 0 and numframes values on axis 1"""
        def fget(self):
            return self._data[:-1]
        def fset(self,x):
            self._data[:-1] = x
        return locals()
    trajectory = property(**trajectory())

    def site_trajectory(self,sites):
        """Returns numpy array with N sites on axis 0 and numframes values on axis 1; indexed by the site label"""
        return self.trajectory[self.site2indx[sites]]

    def site_distribution(self,sites):
        """Returns the distributions for the array <sites>."""
        return self.distribution[self.site2indx[sites]]

    # interface to the super class:
    data = trajectory
    site_data = site_trajectory

class SiteHistogramArray(SiteArray):
    """A collection of histograms (or distributions), one for each site.

    Organized as a big numpy array. The first axis corresponds to the site
    index, the second axis contains the bins. The first and the last bin
    contain counts for data below and above the first (lo_bin) and last bin
    edge (hi_bin).

    :Methods:
    normalize()    create the distribution
    plot()         plot the distribution
    save()         save to a pickled file
    load()         load from a pickled file

    :Attributes:

    Arrays are 0-indexed and correspond to the sites. In order to access data
    with the site label s use SiteAnalysis.site2indx[s] as the index.

    histogram      the histograms as one array of shape (Nsites,Nbins), i.e. histogram[k]
                   is the k-th histogram. NOTE: to get the histogram of site s use
                   histogram[SiteAnalysis.site2indx[s]]
    distribution   normalized histogram, created with normalize()
    normalization  total number of counts for the sites
    edges          bin edges
    midpoints      bin centers
    delta          bin sizes
    site2indx      mapping of site label to an index used in all the arrays
    sites          mapping of index to site label
    """

    def __init__(self,name,sites,maxsites,parameters,properties=None,**kwargs):
        """Initialize the SiteHistogramArray.

        H = SiteHistogramArray(name,sites,maxsites,
                               parameters=dict(Nbins=<int>,lo_bin=<float>,hi_bin=<float>),
                               properties=<dict>)

        :Arguments:
        sites          list of all site labels
        maxsites       highest site label encountered in the trajectory
        parameters     dict of
                         Nbins          number of bins
                         lo_bin         lower edge of first bin
                         hi_bin         upper edge of last bin
        properties     dict to customize the plot (xlabel, legend_fontsize, ...).
        """
        super(SiteHistogramArray,self).__init__(name,sites,maxsites,parameters,
                                                properties=properties,
                                                **kwargs)
        # packaged as a dict for a cleaner interface
        try:
            Nbins,lo_bin,hi_bin = parameters['Nbins'],parameters['lo_bin'],parameters['hi_bin']
        except:
            raise KeyError('parameters must be a dict containing Nbins, lo_bin, hi_bin')

        # +1 histogram for all sites not collected (idx=Nsites),+2 bins (0 and
        # Nbins+1) for outliers. Uses array instead of dict (wastes space but
        # allows direct indexing, eg 'site2indx[sites]'; excluded sites index
        # into the 'overflow' histogram idx=Nsites); unfortunatley, because the
        # site label is used as the index into the array, we need to make the
        # array big enough to encompass ALL site labels that could come up in
        # the trajectory (maxsites). This small memory waste allows us to do fast
        # lookup calculations on whole numpy arrays instead of using slow loops on dicts.

        # All Nsites histograms have the same bins; have to add an extra edge
        # (for digitize), assuming that all bins have the same width
        self.bins = numpy.linspace(lo_bin,hi_bin,Nbins,endpoint=False)  # same as in numpy.histogram
        self.edges = numpy.concatenate((self.bins, [2*self.bins[-1]-self.bins[-2]])) # for digitize()
        self.midpoints = 0.5*(self.edges[1:] + self.edges[:-1])
        self.delta = self.edges[1:] - self.edges[:-1]

        # _histogram holds the raw data. axis 0 = sites, axis 1 = bins
        # i = site2indx[s] --> indexes site s into _histogram[i,:]
        # s = sites[i]      --> index i in histogram corresponds to site s
        self._histogram = numpy.zeros((self.Nsites+1, Nbins+2),dtype=int) # outlier sites&values
        self.distribution = numpy.zeros((self.Nsites, Nbins))             # normalized,no outliers

    def reset(self):
        """Erases all accumulated data."""
        self._histogram *= 0
        self.distribution *= 0

    def digitize(self,site_data):
        """Returns the bin for each value X(site)."""
        # NOTE: site_data[] contains a single value per site
        # (Perhaps should accumulate and digitize larger lists?  E.g. use a
        # buffer array and once one row fills up, accumulate the histograms;
        # this _might_ allow for more efficient histogramming.)
        return [ numpy.digitize([x],self.edges)[0] for x in site_data ]  # too slow? use buffer?

    def add(self,site_data,sites):
        """Add one frame of data to the site histogram.

        site_data     numpy array; axis=0 indexes the sites
                      For scalars, site_data.shape == (N,), with N == len(sites)
                      For vectors, site_data.shape == (N,3).
        sites         numpy array of the site labels that correspond to the slices
                      along axis=0 in site_data
        """
        # For each site in sites, add +1 to the bin corresponding to site_data[site].
        # sites and site_data must have the same shape, i.e. (len(sites),)
        self._histogram[self.site2indx[sites], self.digitize(site_data)] += 1

    def normalize(self):
        """Calculate normalization (including outliers) and distribution."""
        self.normalization = numpy.sum(self._histogram[:-1],axis=1) # includes outliers
        self.normalization[self.normalization == 0] = -1  # fix for normalization of empty histos
        self.distribution = self.histogram/(self.normalization[:,numpy.newaxis] * self.delta)

    def histogram():
        doc = """Numpy array with N sites on axis 0 and Nbin bins on axis 1 (view on _histogram)"""
        def fget(self):
            return self._histogram[:-1,1:-1]
        def fset(self,x):
            self._histogram[:-1,1:-1] = x
        return locals()
    histogram = property(**histogram())

    def site_histogram(self,sites):
        """Returns numpy array with N sites on axis 0 and Nbin bins on axis 1; indexed by the site label"""
        return self.histogram[self.site2indx[sites]]

    def site_distribution(self,sites):
        """Returns the distributions for the array <sites>."""
        return self.distribution[self.site2indx[sites]]

    # interface to the super class:
    data = histogram
    site_data = site_histogram

    def mean(self):
        """Mean for each site <X> (first moment of distribution, <X> = sum x*p(x))."""
        # use self.bins instead of self.midpoints as this is arguably more correct for
        # occupancies (and makes it consistent with occupancy from trajectory)
        return numpy.sum(self.distribution * self.bins, axis=1)

    def std(self):
        """Standard deviation for each site sqrt(<(X-<X>)**2>), from second moment sum x**2*p(x)."""
        EX = self.mean()
        EXX = numpy.sum(self.distribution * self.bins**2, axis=1)
        return numpy.sqrt(EXX - EX*EX)


class SiteAnalysisObservable(object):
    """Base class to implement 'plugin' analysis for site-based trajectory analysis.

    Derive a class, register it in the global observable_registry, and use the
    instance in the compute() method of the Collection class.

    The subclasses override __init__ with the following argument list:

    Observable(collection, histogram_parameters, trajectory=False, **kwargs)

    collection            hop.siteanalysis.Collection to which the observable belongs
                          (uses many attributes for setup)
    histogram_parameters  defines bins and range of histograms
    trajectory            True: also collect time series of the observable, resolved
                                by site i: observable_i(t)
                          False: only compute the histograms (default)
    **kwargs              additional data (subclasses pick what they need)

    The subclass picks anything it needs from the auxiliary kwargs dict and adds it as
    attributes so that its own compute() method can access it; they must silently ignore
    anything they don't need.

    Also override the class-variables
        defaultbinwidth      (for automatic setting of Nbins)
        plot_properties      dict, see below
        name                 same name as used in hop.siteanalysis.ObservablesRegistry (IMPORTANT!!!)

    :Attributes:
    (mainly from <collection>)

    sites        selection of site labels; only on these sites histograms are built
    max_site     the highest site label that occurs in the trajectory
    idcd2ihop    translation between dcd indices to hop indices
                 sites[self.idcd2ihop[[2,3,5]] -> sites of dcd 2,3,5
    data_indices == MDAnalysis.selectAtoms(...).indices() (atom indices of selected atoms)
    data_ts      (pointer to) the data trajectory time step instance
    hop_ts       (pointer to) the hop trajectory time step instance

    (from __init__)
    histogram_parameters = {'Nbin': <number of bins in histogram>, 'lo_bin':<>, 'hi_bin':<>}
    plot_properties      = {'xlabel': <>, 'title': <>, ...}
    """

    defaultbinwidth = 0.5    # set default bindwidth for each subclass
    name = None              # name used in ObservablesRegistry

    def __init__(self,collection,histogram_parameters,plot_properties=None,trajectory=False,**kwargs):
        """Set up an observable (a class describing how a histogrammed site property is computed).

        o = SiteAnalysisObservable(collection,histogram_parameters,plot_properties)

        :Arguments:
        collection            hop.siteanalysis.Collection instance
        histogram_parameters  {'Nbin': <number of bins in histogram>, 'lo_bin':<float>,
                               'hi_bin':<float>}
        trajectory            False: do not set up SiteTrajectoryArray (i.e. do not store the
                              whole site trajectory) or True (allocate enough space to hold
                              the _whole_ trajectory!)
        plot_properties       {'xlabel': <string>, 'title': <string>, ...}
        """
        super(SiteAnalysisObservable,self).__init__()  # necessary for classes with super()
        self.collection = collection
        self.sites = collection.siteselection
        self.all_sites = collection._all_sites       # all 'default' sites in the trajectory (from density)
        self.max_site = collection._all_sites.max()  # largest site label
        self.idcd2ihop = collection.idcd2ihop
        self.data_indices = collection.data_indices
        self.data_ts = collection.data_ts        # pointers to changing Timestep
        self.hop_ts = collection.hop_ts

        if plot_properties is None:
            plot_properties = self.plot_properties  # self.plot_properties should be set in subclass

        # set up site histograms
        if histogram_parameters['hi_bin'] <= histogram_parameters['lo_bin']:
            raise ValueError("histogram parameters must be hi_bin > lo_bin, but "
                             "lo_bin=%(lo_bin)g and hi_bin=%(hi_bin)g." % histogram_parameters)
        if not histogram_parameters['Nbins']:    # fill in default for Nbins
            histogram_parameters['Nbins'] = int(\
                (histogram_parameters['hi_bin'] - histogram_parameters['lo_bin'])/self.defaultbinwidth)

        self.histograms = SiteHistogramArray(self.name,self.sites,self.max_site,
                                             histogram_parameters,
                                             properties=plot_properties)

        # set up site trajectory arrays (if requested)
        self.has_trajectories = trajectory
        if trajectory:
            assert collection.data_dcd.numframes == collection.hop_dcd.numframes
            trajectory_parameters = {'numframes': collection.hop_dcd.numframes,
                                     'totaltime': trajectory.totaltime(collection.hop_dcd,'ps'),
                                     'delta_t':   trajectory.delta_t(collection.hop_dcd,'ps')}
            self.trajectories = SiteTrajectoryArray(self.name,self.sites,
                                                    self.max_site,trajectory_parameters,
                                                    properties=plot_properties)

    def compute(self,**kwargs):
        """Compute the histograms (and trajectories) for selected sites for the given time step."""
        self.build_histograms(**kwargs)
        if self.has_trajectories:
            self.build_trajectories(**kwargs)  # may depend on values computed in build_histograms()

    def build_histograms(self,**kwargs):
        # see Orbit.build_histograms() with self._compute_distance() as an example
        """Override this method. It must silently ignore any keyword args that it does not use."""
        raise NotImplementedError("Internal Error: The build_histograms() method of class %r was not implemented." %
                                  self.__class__)

    def build_trajectories(self,**kwargs):
        """Override this method. It must silently ignore any keyword args that it does not use."""
        raise NotImplementedError("Internal Error: The build_trajectorie() method of class %r was not implemented." %
                                  self.__class__)

    def reset(self):
        """Clear all computed data."""
        self.reset_histograms()
        self.reset_trajectories()

    def reset_histograms(self):
        try: self.histograms.reset()
        except AttributeError: pass

    def reset_trajectories(self):
        try: self.trajectories.reset()
        except AttributeError: pass

    # utility methods:
    # Try to avoid reallocating numpy arrays: use assignments to buffer arrays
    # self._NAME (which are allocated in 'try: ... except AttributeError: ...' constructs)

    def _get_sites(self,hop_pos):
        """Make selected site labels of the time step available in self._sites[]. They
        are arranged in the same order as the SiteHistogramArray and the
        SiteTrajectoryArray."""
        # (Perhaps this should be made a class variable (cache) and only computed once
        # for ALL observables?)
        try:
            # site label for each selected data atom
            self._sites[:] = hop_pos.astype(int)[self.idcd2ihop]
        except AttributeError:  # allocate for the first time
            self._sites = hop_pos.astype(int)[self.idcd2ihop]

class Distance(SiteAnalysisObservable):
    """Observable: distance of atom from the centre of a site while on the site.

    o = Distance(Collection,histogram_parameters,site_centers=site_centers)
    """

    defaultbinwidth = 0.1   # 0.1 A because distance is typically 0 <= d < 4 A
    plot_properties = {'xlabel':"distance [Angstrom]",'title':"Distance from site"}
    name = "distance"

    def __init__(self,*args,**kwargs):
        super(Distance,self).__init__(*args,**kwargs)
        self._init_site_centers(kwargs['density'])

    def build_histograms(self,**kwargs):
        self._compute_distance(self.hop_ts._x)  # use instantaneous site s(t)

    def build_trajectories(self,**kwargs):
        # TODO: Must reduce data for site, eg average, max, min, median, ...
        # THIS IS WRONG (just chooses one of the values)
        raise NotImplementedError('Current implementation of distance site traj is broken.')
        self.trajectories.add(self._dist,self._sites)

    # auxiliary functions that can be directly used by subclasses
    # (_compute_distance() shows an example implementation of build_histograms() )
    def _compute_distance(self,hop_pos):
        """Accumulate site histograms for distance site--coord.

        pos          sites from the hopping trajectory

        Uses self.site_centers centers of all sites (i.e. site_properties.center))
        which must be initialized with _init_site_centers(density).
        """
        # allocate arrays when needed (but only once so that we can use fast assignment)
        self._get_sites(hop_pos)
        try:
            try:
                self._distvec[:] = self.site_centers[self._sites] - self.data_ts._pos[self.data_indices]
            except AttributeError:
                self._distvec = self.site_centers[self._sites] - self.data_ts._pos[self.data_indices]
        except IndexError,err:
            raise IndexError("Encountered a site in the hop trajectory that was not listed "
                             "in the density ("+str(err)+"). Check the consistency of the "
                             "input data.\n")
        try:
            self._dist[:] = numpy.sqrt(numpy.sum(self.distvec*self.distvec,axis=1)) # observable
        except AttributeError:
            self._dist = numpy.sqrt(numpy.sum(self._distvec*self._distvec,axis=1)) # observable
        self.histograms.add(self._dist,self._sites)           # accumulate in histograms for each site

    def _init_site_centers(self,density):
        # site centers: reformat array:
        # (1) 2D array, (2) [NaN,NaN,NaN] instead of None
        site_centers = density.site_properties.center.copy()
        site_centers[SITELABEL['interstitial']] = numpy.array(3*[numpy.NaN])
        self.site_centers = numpy.array(site_centers.tolist())  # shape == (num_sites,3)

class Orbit(Distance):
    """Observable: distance of atom from the centre of the site while orbiting the site.

    o = Orbit(Collection,histogram_parameters,site_centers=site_centers)

    (A particle orbits a site if it has been on the site at a previous time step but has
    not entered another site since.)
    """

    defaultbinwidth = 0.2    # A, typical ranges for orbit 0 <= o < 15 A
    plot_properties = {'xlabel':"distance [Angstrom]",
                       'title':"Site orbit",
                       }
    name = "orbit"

    def build_histograms(self,**kwargs):
        self._compute_distance(self.hop_ts._y)  # use orbit o(t) from hoptraj

class Occupancy(SiteAnalysisObservable):
    """Observable: occupancy (number of particles on the site).

    o = Occupancy(Collection,histogram_parameters,trajectory=False)

    histogram_parameters    dict with lo_bin, hi_bin, [Nbins]
    trajectory              True: also extract the site occupancy as a function of
                                  time (or frame), s_i(t)
                            False: only do the histogram
    """

    defaultbinwidth = 1   # one bin for each particle per site (only integers)
    #plot_properties = {'xlabel':r"$\nu_\mathrm{site}$",'title':"Site occupancy distribution"}
    plot_properties = {'xlabel':r"$\nu$",
                       'title':"Site occupancy distribution",
                       'xticks_indexed':True}
    name = 'occupancy'

    def __init__(self,*args,**kwargs):
        super(Occupancy,self).__init__(*args,**kwargs)
        # quick and dirty: histogram on ALL sites
        # Add one additional bin at the higher end which holds the outliers at the upper
        # end; the ones at the lower end are automatically discarded
        hmin, hmax = self.all_sites.min(), self.all_sites.max()+1
        self._bins = numpy.linspace(hmin,hmax,num=hmax-hmin+1,endpoint=True)
        # histogrammed_sites[]: index in bins --> site label (includes unselected sites!)
        self.histogrammed_sites = self._bins.astype(int)[:-1]  # remove last bin = outliers
        self._histo = self._occupancy_histogram([0])  # initialize histogram
        self._histo *= 0

    def _occupancy_histogram(self,data):
        """Return the histogram with the outlier bin discarded."""
        return numpy.histogram(data,bins=self._bins)[0][:-1]

    def _compute_occupancy(self,hop_pos,**kwargs):
        self._get_sites(hop_pos)   # --> self._sites
        # X(s) where s: selected site, X: number of particles on site
        self._histo[:] = self._occupancy_histogram(self._sites)
        self.histograms.add(self._histo, self.histogrammed_sites)    # accumulate in histograms for each site

    def build_histograms(self,**kwargs):
        self._compute_occupancy(self.hop_ts._x)  # x: occ, y: orbit_occup

    def build_trajectories(self,**kwargs):
        # Make sure that this runs AFTER build_histograms(): needs _histo!
        self.trajectories.add(self._histo, self.histogrammed_sites)

class OrbitOccupancy(Occupancy):
    """Observable: occupancy, based on orbit.

    o = OrbitOccupancy(Collection,histogram_parameters)
    """

    #plot_properties = {'xlabel':r"$\nu_\mathrm{site}$",'title':"Site orbit occupancy distribution"}
    plot_properties = {'xlabel':r"$\nu$",'title':"Site orbit occupancy distribution"}
    name = "orbitoccupancy"

    def build_histograms(self,**kwargs):
        self._compute_occupancy(self.hop_ts._y)  # x: occ, y: orbit_occup


# dict of SiteAnalysisObservable classes
ObservablesRegistry = {'distance':Distance,
                       'orbit':Orbit,
                       'occupancy':Occupancy,
                       'orbitoccupancy':OrbitOccupancy,
                       }

class Collection(utilities.Saveable):
    """A Collection of observables to be calculated for a hopping trajectory.
    (the central object of the siteanalysis module)

      C = Collection(...)
      C.add_observable(...)
      C.add_observable(...)
      ...
      C.compute()
      C.plot()

    The general outline is:

    (1) Load dcd with selection that contains the atoms in the hopping
    trajectory. Note that each residue in this selection must have exactly one
    atom in the hopping trajectory (or the selection must be identical to the
    hopping selection); also load the corresponding hopping trajectory:

      C = Collection(selection=universe.selectAtoms('name OH2'),
                     hoptraj=HoppingTrajectory(filename='hoptraj'),
                     density=Density(filename='water.pickle'))

    (2) Register observables with the Collection:

      C.add_observable('distance', lo_bin=0, hi_bin=3)

    Note that only observables listed in ObservablesRegistry can be used. In
    order to use new observables derive from the class SiteAnalysisObservable
    and look at the source code for the Distance class as an example.

    (3) At each time step, calculate the observable X for all atoms a in the
    selection and associate the value X(a) with the site s(a).

      C.compute()

    :Methods:
    add_observable()   Add another observable to the Colletion

    compute()          Do the analysis. Specify the observables here.
    imshow()           Plot 2D histograms.
    plot()             Plot histograms.
    load()             Instantiate from pickled file.
    save()             Pickle histograms to file.
    """

    _saved_attributes = ['siteselection','site_properties','_all_sites',
                         'SiteHistograms','SiteTrajectories',
                         ]

    def __init__(self,selection=None,hoptraj=None,density=None,
                 include_sites='default',exclude_sites=['default','bulk'],
                 filename=None, verbosity=None):
        """Set up analysis on a per-site basis.

        C = Collection(selection=universe.selectAtoms('name OH2'),
                  hoptraj=HoppingTrajectory(filename='hoptraj'),
                  density=Density(filename='water.pickle'),
                  include_sites='default',exclude_sites=['default','bulk'],
                  filename=None, verbosity=None)

        :Arguments:
        selection          MDAnalysis.selectAtoms selection (from psf+dcd)
        hoptraj            hop.trajectory.HoppingTrajectory instance (from psf+dcd)
        density            hop.sitemap.Density instance with a sitemap defined
        include_sites      any valid argument to
        exclude_sites      Density.site_labels(include=include_sites,exclude=exclude_sites)
                           The default accumulates data for all the sites in the trajectory
                           except for the bulk site.
                           Histograms will only be produced for those sites. Note that
                           the calculation will not be much faster for fewer sites, only memory
                           requirements for the histograms will be lower.
        filename           filename to write to or read histograms from
        """
        set_verbosity(verbosity)
        if not (selection is None and hoptraj is None and density is None):
            # calculate anew from raw data
            self.universe = selection.atoms[0].universe
            self.selection = selection
            self.data_dcd = self.universe.trajectory
            self.hop_dcd = hoptraj.hoptraj
            self.data_ts = self.data_dcd.ts   # pointers to the changing time step object
            self.hop_ts = self.hop_dcd.ts
            self.density = density  # needed for specific observables AND for ALL the sitelabels
            self.siteselection = density.site_labels(include=include_sites,exclude=exclude_sites)
            self.site_properties = density.site_properties
            self._all_sites = density.site_labels('default') # all the proper sites in a trajectory
            self.site_labels_cache = {'all_sites': density.site_labels('default'),
                                      'subsites': density.site_labels('subsites'),
                                      }
            # indices in ts._pos of the selected atoms (used for 'ts._pos[data_indices]' which is
            # slightly faster than 'self.selection.coordinates()' because it avoids a copy)
            self.data_indices = self.selection.indices()  # needed for Observables compute()
            self.SiteHistograms = {} #IntrospectiveDict()     # will hold the results after compute()
            self.SiteTrajectories = {} #IntrospectiveDict()
            self.observable_plugins = {} #IntrospectiveDict() # observable (name:instance) for compute()
            if filename is not None:
                self.filename(filename,set_default=True)

            # sanity checks
            if self.data_dcd.numframes <> self.hop_dcd.numframes:
                raise ValueError("Data trajectory and hopping trajectory differ in length: "
                                 "%d vs %d frames" % (self.data_dcd.numframes,self.hop_dcd.numframes))

            # compute index arrays between sites, atoms, and histograms
            hopSR2A = {}
            for ihop,atom in enumerate(hoptraj.group.atoms):
                sr_index = (atom.segment.name,atom.residue.id)
                if sr_index in hopSR2A:
                    raise ValueError("Each residue should only have a single atom in the "
                                     "HoppingTrajectory. Detected a collision:\n"
                                     "Previous atom index " + str(hopSR2A[sr_index]) + ":"\
                                         +str(hoptraj.group[hopSR2A[sr_index]])+"\n"
                                     "New atom index " + str(ihop) + ":    "\
                                         +str(atom)+"\n")
                hopSR2A[sr_index] = ihop
            idcd2ihop = []  # map index in data dcd to index in hop (many-to-one!)
            for idcd,atom in enumerate(selection.atoms):
                sr_index = (atom.segment.name,atom.residue.id)
                try:
                    idcd2ihop.append(hopSR2A[sr_index])
                except KeyError:
                    raise KeyError("Atom with index "+str(idcd)+" "+str(atom)+\
                                       " in the data trajectory belongs to a "
                                   "residue that has no atom in the hop trajectory. Change the "
                                   "selection appropriately.")
            self.idcd2ihop = numpy.array(idcd2ihop)
            del idcd2ihop
            del hopSR2A
        elif filename is not None:
            self.load(filename)
        else:
            raise ValueError("Not enough data in argument list. Either supply selection, "
                             "hoptraj, and density or load from a pickled file <filename>.pickle")

    def add_observable(self,name,lo_bin,hi_bin,Nbins=None,**kwargs):
        """Add observable <name> to the Collection.

        add_observable(<name>,<lo_bin>,<hi_bin>,Nbins=<Nbins>,**kwargs)

        The observable <name> is added to the collection so that <name> is
        histogrammed for all selected sites over the trajectory when the
        Collection.compute() method is run. The number of bins, <Nbins>, and the
        lower edge of the first bin, <lo_bin>, and the upper edge of the last
        bin, <hi_bin>, must be supplied.

        :Arguments:
        name        name of observable (see hop.siteanalysis.ObservablesRegistry)
        lo_bin      value for the lower edge of the first bin
        hi_bin      value for the upper edge of the last bin
        Nbins       number of bins (None sets it to a default value, depending on <name>)
        kwargs      keyword arguments specifying additional parameters for the particular
                    observable (see docs for the observables)
        """
        histogram_parameters = {'Nbins':Nbins,'lo_bin':lo_bin,'hi_bin':hi_bin}
        try:
            self.observable_plugins[name] = \
                ObservablesRegistry[name](self,histogram_parameters,density=self.density,**kwargs)
        except KeyError:
            raise ValueError('Observable '+str(name)+' not implemented. Choose one of '
                             +str(ObservablesRegistry.keys()))

    def del_observable(self,names=None):
        """Remove observable plugins <name> or [<name1>,<name2>,...] from the list of computed observables."""
        for o in self._observables(names):
            del self.observable_plugins[name]

    def reset_observable(self,names=None):
        """Reset observable histograms to initial empty values; None means 'all'."""
        msg(5,'Resetting all observables so that we can accumulate new data.\n')
        for o in self._observables(names):
            o.reset()

    def _observables(self,names=None):
        try:
            if names is None:
                return self.observable_plugins.values()
            elif iterable(names):
                return [self.observable_plugins[name] for name in names]
            else:
                return [self.observable_plugins[names]]
        except KeyError:
            raise ValueError("Only observables defined in %r can be used as <observables>." %
                             self.observable_plugins.keys())

    def compute(self,observables=None,start=None,stop=None,step=1,progress_interval=10):
        """Compute histograms of the observables.

        :Arguments:
        observables   list of names of observables to compute; None means 'all'
        start         0...numframes-1  NOTE: 0-based frame indexing
        stop          1...numframes
        step          1
        progress_interval   10 (print progress every <interval> steps)

        :Attributes:
        The method populates the dict SiteHistograms with the results.

        SiteHistograms Dict of the histograms, one for each observable. Each
                       entry o contains one SiteHistogramArray o, which holds
                       an array of shape (Nsites,Nbins), i.e. o.histogram[k] is
                       the k-th histogram. NOTE: to get the histogram of site s
                       use o.site_histogram(s)
        SiteTrajectories
                       If the observable has the 'trajectory=True' flag enabled, the
                       value of the observable for each site and time frame is stored.
        """
        # sanity check
        try:
            self.data_dcd.numframes
            self.hop_dcd.numframes
            self.density.site_properties.center
        except:  # typically triggered when Collection was recreated from a pickled file
            raise MissingDataError("Cannot compute() histograms because trajectory "
                                       "data is missing.")
        if start is None: start = 0
        elif start < 0:   start = self.data_dcd.numframes + start
        if stop is None:  stop = self.data_dcd.numframes
        elif stop < 0:    stop = self.data_dcd.numframes + stop + 1
        numframes = (stop - start + 1)/step

        _observables = self._observables(observables)
        self.reset_observable(observables)
        msg(3,'Computing observables %r\n' % [o.name for o in _observables])

        # loop through all frames and add data (compute) for each frame
        for iframe,percentage in self._synchronous_trajectories_iter(start,stop,step):
            for o in _observables:
                o.compute()
            if iframe % progress_interval == 0:
                msg(3,"Data/Hop frame %4d  [%4.1f%%]\r" % (iframe, percentage))
        msg(3,"Data/Hop frame %4d  [%4.1f%%]\n" % (iframe, percentage))
        msg(5,'Creating normalized distributions.\n')
        for o in _observables:
            o.histograms.normalize()
        # only keep SiteHistogramArrays and SiteTrajectoryArrays
        self.SiteHistograms.update(dict([(o.name,o.histograms) for o in _observables]))
        self.SiteTrajectories.update(dict([(o.name,o.trajectories) for o in _observables
                                           if o.has_trajectories]))


    def plot(self,observable, *args,**kwargs):
        """Plot all histograms or trajectories (or a selection of sites).

        plot(observable, filename=None,sites=None,format='pdf',with_legend=True)

        :Arguments:
        observable     name of the observable to plot
        filename       output graphic filename (extension determines format or use format=...)
        sites          None: everything in the Collection; or a list of numeric site labels
        trajectory     False: plot histogram, True: plot trajectory (if computed)
        with_legend    True or False
        """
        self._plot_with('plot',observable,*args,**kwargs)

    def imshow(self, observable, *args,**kwargs):
        """Plot all histograms or trajectories (or a selection of sites) as a 2D density map.

        imshow(observable, filename=None,sites=None,format='pdf')

        :Arguments:
        observable     name of the observable to plot
        filename       output graphic filename (extension determines format or use format=...)
        sites          None: everything in the Collection; or a list of numeric site labels
        trajectory     False: plot histogram, True: plot trajectory (if computed)
        normalize      True: distribution; False: histogram
        format         determine output format if filename is derived from the
                       default filename
        **kwargs       keyword arguments are passed directly to pylab.imshow()
        """
        self._plot_with('imshow',observable,*args,**kwargs)

    def _plot_with(self,methodname,observable,*args,**kwargs):
        trajectory = kwargs.pop('trajectory',False)
        if trajectory:
            results = self.SiteTrajectories
        else:
            results = self.SiteHistograms
        try:
            sitearray = results[observable]
        except KeyError:
            raise ValueError("'observable' must be one of "+str(self.results.keys()))
        plotfunctions = {'plot':sitearray.plot,
                         'imshow':sitearray.imshow,
                         }
        try: filename = args[0]
        except:
            try:  filename = kwargs['filename']
            except: # no filename in argument list: use our default
                kwargs['filename'] = self.filename()

        plotfunctions[methodname](*args,**kwargs)

    def _next_frame(self,start=None):
        if start is None:
            data_ts = self.data_dcd.next()
            hop_ts = self.hop_dcd.next()
        else: # or should this be coded with repeated next() to avoid random access bugs?
            data_ts = self.data_dcd[start]
            hop_ts = self.hop_dcd[start]
        assert(data_ts.frame == hop_ts.frame)

    def _synchronous_trajectories_iter(self,start,stop,step=1):
        """Go through hop and data frames by only using the next() method thus
        avoiding possible problems with random-access of big trajectories or non-RA
        formats such as XTC."""
        total_frames = stop - start
        self._next_frame(start)          # advance to first frame
        for iframe in xrange(start,stop,step):
            try:
                yield iframe,100.0*(iframe-start+1)/total_frames
                [self._next_frame() for j in xrange(step)]   # advance at the end of loop!
            except IOError:
                raise StopIteration

    def interactive(self):
        """Turn SiteHistograms and SiteTrajectories into introspection-friendly objects cls.H and cls.T
        """
        self.H = IntrospectiveDict(self.SiteHistograms)
        self.T = IntrospectiveDict(self.SiteTrajectories)
        msg(0,"H contains %r\n" % self.H.keys())
        msg(0,"T contains %r\n" % self.T.keys())

#
#----------------------------------------------------------------------
#

class Test(object):
    def __init__(self,psf='inp/ifabp_water.psf',dcd='trj/rmsfit_ifabp_water_1.dcd',
                 density="analysis/water",hoptrj="trj/hoptrj",with_observables=True):
        """Test typical setup for SiteAnalysis. Supply filenames for
        psf and dcd (the data trajectory), density (a pickled Density
        object with a sitemap), and the name of the hopping trajectory
        (it assumes that the hopping psf has the same name with the
        psf prefix)"""
        import MDAnalysis
        #utilities.matplotlib_interactive(False)
        import hop.sitemap

        self.u = MDAnalysis.Universe(psf,dcd)
        self.wat = self.u.selectAtoms('name OH2')
        self.hoptraj = trajectory.HoppingTrajectory(filename=hoptrj)
        self.dens = hop.sitemap.Density(filename=density)
        self.A = Collection(self.wat,self.hoptraj,self.dens)
        if with_observables:
            self.A.add_observable('distance', lo_bin=0, hi_bin=3)
            self.A.add_observable('orbit', lo_bin=0, hi_bin=15)
            self.A.add_observable('occupancy', lo_bin=0, hi_bin=4, trajectory=True)
            self.A.add_observable('orbitoccupancy', lo_bin=0, hi_bin=10, trajectory=True)

        print "defined observables: %r" % self.A.observable_plugins.keys()
        print "data acquistion for observables: T.A.compute(...)"
        print "plot results: T.A.plot(T.A.sites)"
        print "ATTENTION: is the density really correct? Should be 'remapped' for holo"\
            " density name = %(density)s" % locals()

        # self.A.compute()
        # self.A.plot(self.A.sites)

def _fix_old(filename='analysis/siteanalysis'):
    A = hop.siteanalysis.Collection(filename=filename)
    for n,x in A.SiteHistograms.items():
        p = x.default_plotproperties.copy()
        p.update(x.properties)
        x.properties = p
        print n
        print x.properties

    A.save()

    A.imshow('distance',vmax=None,filename='figs/site_distance2.pdf')
    A.imshow('orbit',vmax=None,filename='figs/site_orbit2.pdf')
    A.imshow('occupancy',filename='figs/site_occupancy.pdf')
    A.imshow('orbitoccupancy',filename='figs/site_orbitoccupancy.pdf')
