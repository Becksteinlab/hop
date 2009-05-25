# $Id$

"""A collection of functions and classes to extract statistics and plot
histograms. Use this as examples how to write your own.
"""

import hop.graph, hop.sitemap, hop.utilities
from hop.utilities import msg,set_verbosity,verbosity
import numpy
import pylab
import cPickle
import pprint
import os, os.path

from hop.utilities import sorted


pp = pprint.PrettyPrinter(indent=4)


class DensityAnalysis:
    def __init__(self,densities,reference,dir="./figs",
                 bulkdensities=None,refbulkdensity=None,
                 bulkname='bulk',verbosity=3):
        """Analyze densities.

        :Arguments:
        densities      list of densities (instances or filenames)
        reference      reference density (instance or filename)
        bulkdensities  needed to rebuild densities with scan; same
                       sequence as densities
        refbulkdensity ...
        bulkname       eg 'bulk': searches for '<bulkname>.pickle' in
                       the same directory as the density
        dir            basedir for graphs
        verbosity=3    chattiness

        :Methods:
        show()         basic stats
        scan()         scan statistics for a range of rho_cut
        save()         pickle scan results
        load()         load scan results
        """
        set_verbosity(verbosity)

        self.densities = hop.utilities.easy_load(densities,hop.sitemap.Density,'stats')
        self.reference = hop.utilities.easy_load(reference,hop.sitemap.Density,'stats')
        try:
            self.n = len(self.densities)   # densities are numbered 0 .. n-1
        except TypeError:
            self.densities = [self.densities]
            self.n = 1

        if bulkdensities and refbulkdensity:
            self.bulkdensities = hop.utilities.easy_load(bulkdensities,hop.sitemap.Density,'stats')
            self.refbulkdensity = hop.utilities.easy_load(refbulkdensity,hop.sitemap.Density,'stats')
        elif bulkname:
            self.refbulkdensity = hop.utilities.easy_load(
                os.path.join(os.path.dirname(self.reference.filename()),bulkname),
                hop.sitemap.Density,'stats')  
            self.bulkdensities = [ hop.utilities.easy_load(
                    os.path.join(os.path.dirname(d.filename()),bulkname),
                    hop.sitemap.Density,'stats') for d in self.densities ]
        else:
            self.bulkdensities = None
            self.refbulkdensity = None

        if self.bulkdensities and type(self.bulkdensities) is not list:
            self.bulkdensities = [self.bulkdensities]


        # D contains raw data, indexed by density or 'reference'
        self.D = dict( [(k,{}) for k in xrange(self.n)] )
        self.D['reference'] = {}
        # S contains all the statistics
        self.S = dict( [(k,x.stats(data=self.D[k])) for k,x in enumerate(self.densities)] )
        self.S['reference'] = self.reference.stats(data=self.D['reference'])
        # scanner is used for scan()
        self.scanner = DensityScanner(self)
        # default pdf file names
        self.pdf = {'scan':'scan.pdf',
                    }
        for k,v in self.pdf.items():
            self.pdf[k] = os.path.join(dir,v)

    filename = hop.utilities.filename_function

    def all(self):
        self.show()

    def show(self):
        print "Statistics S"
        print "------------"
        pp.pprint(self.S)

        print "Data in D"
        print "---------"
        pp.pprint(sorted(self.D.keys()))

    def scan(self,cutoffs,merge=True):
        """Recompute density stats for various cutoffs  (thresholds)
        in self.scanstats.

          scan([1.5,2.0, 2.72, 3.0],merge=True)

        merge        True: add/replace to scanstats, False: erase
        See also DensityAnalysis.scanner.scanarrays for a dict of
        numpy.recarrays of the data.
        """
        self.scanner.scan(cutoffs,merge)

    def save(self,filename):
        """Pickle the results from scan()."""
        self.scanner.save(filename)
    def load(self,filename):
        """Load the pickled results from a scan."""
        self.scanner.load(filename)
    def plotscan(self,filename=None,properties=None,fignumber=1):
        """Plot results of a scan(). 
        
        plotscan('figs/scan.pdf')

        See scanner.plot for docs."""
        if filename is None:
            filename = self.pdf['scan']
        self.scanner.plot(filename,properties=properties,fignumber=fignumber)        

class DensityScanner(hop.utilities.Saveable):
    # use scanstats for archival because it is more flexible than the record and can be easily updated
    _saved_attributes = ['scanstats']
    def __init__(self,densityAnalysis,with_densities=True):
        self.scanstats = {}         # results from scan()
        self.with_densities = with_densities
        if with_densities:
            self.reference = densityAnalysis.reference
            self.refbulkdensity = densityAnalysis.refbulkdensity
            self.densities = densityAnalysis.densities
            self.bulkdensities = densityAnalysis.bulkdensities

    def scan(self,cutoffs,merge=True):
        if not self.with_densities or not (self.bulkdensities and self.refbulkdensity):
            # need bulk densities because redrawing the map at a new
            # cutoff completely erases the manually inserted bulk
            # density
            raise ValueError("Densities (including bulk) are required for a scan.")
        if not merge:
            self.scanstats = {}    # indexed by cutoff, list of stats
        
        self._check_and_remap([self.refbulkdensity])
        self._check_and_remap(self.densities)
        self._check_and_remap(self.bulkdensities)

        try:
            iter(cutoffs)
        except TypeError:
            cutoffs = [cutoffs]
        for rhocut in cutoffs:
            msg(3, "rhocut = %5.2f: " % rhocut)
            r = self.reference
            r.map_sites(rhocut)
            r.site_insert_bulk(self.refbulkdensity)
            for k,d in enumerate(self.densities):
                msg(3,"%d/ref " % k)
                d.map_sites(rhocut)
                d.site_insert_bulk(self.bulkdensities[k])
                d.find_equivalence_sites_with(r)
                self.scanstats[rhocut] = {'reference':r.stats(),k:d.stats()}
            msg(3,"\n")
        self._scan2recarray()
        msg(3, "Scan complete. See scanstats and scanarrays attributes.")

    def _check_and_remap(self,densityList):
        """check that all densities are compatible with the reference
        density (densities are NOT saved)"""
        for k,dens in enumerate(densityList):
            if dens.map.shape != self.reference.map.shape:
                msg(1, "Must remap density "+str(dens)+" to the reference "+str(self.reference)+"\n")
                dens_remapped = hop.sitemap.remap_density(dens,self.reference) # copy
                densityList[k] = dens_remapped    # assignment (no idea why it has to be that clumsy)  

    def load(self,fn,merge=True):
        super(DensityScanner,self).load(fn,merge=merge)
        self._scan2recarray()

    def _scan2recarray(self):
        # unwrapping/assembling table from stats
        stats = self.scanstats
        xkeys = sorted( stats.keys() )
        ids = sorted( stats[xkeys[0]].keys() )
        fieldnames = sorted( stats[xkeys[0]][ids[0]].keys() )
        self.scanarrays = {}
        for id in ids:
            a = {}
            for y in fieldnames:
                a[y] = [stats[x][id][y] for x in xkeys]
            self.scanarrays[id] = numpy.rec.fromarrays(
                [a[y] for y in fieldnames],names=fieldnames)

    __plot_properties_default = {'sites': {'ylim': None,
                                           'ylabel':'',
                                           },
                                 'volume': {'ylim': (0,10),
                                            'ylabel':'site volume/Angstrom^3',
                                            },
                                 'occupancy': {'ylim': (0,0.8),
                                               'ylabel': 'occupancy from rho',
                                               'xlabel':'rho_cutoff',
                                               },             
                                 }

    def plot(self,fn=None,idens=0,functions='all',properties=None,fignumber=1):
        """Plot density statistics against rho_cut for reference (black) and density 0 (red).

        plot(filename,properties=<dict of dicts>)

        Plot various functions of the density cut-off rho_cut. Current
        functions are 'sites', 'volume', 'occupancy', or 'all'.

        Plots can be customized by using the properties dict. To
        change the ylim and add an title to the sites graph, use

           properties = {'sites': {ylim':(0,220),'title':'number of sites'}}

        :Arguments:
        fn            file name for the file; default is scan.pdf. Suffix determines file type.
        idens         number of density plot; the first one is 0 in self.scanarrays[].
        functions     list of function names or 'all'
        properties    dict1 of dicts; keys1: sites, volume, occupancy;
                      keys2: any matplotlib settable property, values2: appropriate values
        fignumber     pylab figure number
        """
        import pylab
        available_functions = ['sites', 'volume', 'occupancy']
        plot_functions = []
        if functions is 'all':
            plot_functions = available_functions
        else:
            for f in functions:
                if f in available_functions:
                    plot_functions.append(f)
                else:
                    raise ValueError('No function '+str(f)+' is available, only '+
                                     str(available_functions))
        props = self.__plot_properties_default.copy()
        if type(properties) is dict:
            for k,userparams in properties.items():
                try:
                    props[k].update(userparams)
                except KeyError:
                    props[k] = userparams

        r = self.scanarrays['reference']
        d = self.scanarrays[idens]
                
        pylab.figure(fignumber)        
        pylab.clf()
        pylab.axes(axisbg='w',frameon=False)

        par = props['sites']
        pylab.subplot(221)
        pylab.plot(r.rho_cut,r.N_equivalence_sites,'ko',label="N_equiv")
        pylab.plot(r.rho_cut,r.N_sites,'kx',label="N_sites")
        pylab.plot(r.rho_cut,d.N_sites,'rx')
        pylab.ylim(par['ylim'])
        #pylab.legend()

        # refine err bar plots
        alpha = 1  # acts on the line not on the errorbar markers
        capsize=0

        par = props['volume']
        pylab.subplot(222)
        pylab.errorbar(r.rho_cut,r.site_volume_avg,yerr=r.site_volume_std,
                       color='k',capsize=capsize,alpha=alpha,label="volume")
        pylab.errorbar(r.rho_cut,d.site_volume_avg,yerr=d.site_volume_std,
                       color='r',capsize=capsize,alpha=alpha,label="volume")
        pylab.ylim(par['ylim'])
        pylab.ylabel(par['ylabel'])
        # [pylab.set(pylab.gca(), k, v) for k,v in par.items()]

        par = props['occupancy']
        pylab.subplot(223)
        pylab.errorbar(r.rho_cut,r.site_occupancy_rho_avg,yerr=r.site_occupancy_rho_std,
                       color='k',capsize=capsize,alpha=alpha)
        pylab.errorbar(r.rho_cut,d.site_occupancy_rho_avg,yerr=d.site_occupancy_rho_std,
                       color='r',capsize=capsize,alpha=alpha)
        #pylab.legend()
        pylab.ylim(par['ylim'])
        pylab.xlabel(par['xlabel'])
        pylab.ylabel(par['ylabel'])

        pylab.savefig(fn)
        msg(1, "Graph saved to %s\n" % fn)

        pylab.close(fignumber)

class HopgraphAnalysis:
    """Comprehensive analysis of an annotated hop graph."""
    def __init__(self,hopgraph,dir=".",verbosity=3):
        """Analyse hopgraph.
        
        a = HopgraphAnalysis(hopgraph)

        The show() method prints statistics on the HoppingGraph and histograms()
        produces a number of plots as pdf files in the current directory.

        :Arguments:
        hopgraph       can be the name of a pickled file or a HoppingGraph
                       instance
        dir            save figures in this directory
        verbosity=3    chattiness

        :Attributes:
        S              statistics dictionary (see keys for explanation)
        D              raw data dictionary
        
        :Methods:
        all()          show() and histograms()
        show()         print stats
        histograms()   produce histograms
        """
        set_verbosity(verbosity)
        h = hop.utilities.easy_load(hopgraph,hop.graph.HoppingGraph,'stats')
        self.pdf = {'degree':'degree.pdf',
                    'occupancy':'occupancy.pdf',
                    'lifetime':'lifetime.pdf',
                    }        
        self.D = {}
        self.S = h.stats(data=self.D)    # pulls out all the statistics from graph
        for k,v in self.pdf.items():
            self.pdf[k] = os.path.join(dir,v)

    def all(self):
        self.show()
        self.histograms()

    def show(self):
        print "Statistics S"
        print "------------"
        pp.pprint(self.S)

        print "Data in D"
        print "---------"
        pp.pprint(sorted(self.D.keys()))

    def histograms(self):
        self._degree_histo()
        self._lifetime_histo()
        self._occupancy_histo()

    def _degree_histo(self,numfig=1):
        def degreehisto(a):
            return  numpy.histogram(a,bins=numpy.arange(0,numpy.max(a)))

        msg(1,"figure %(degree)s: node degree histograms " % self.pdf)
        pylab.figure(numfig)
        pylab.clf()
        pylab.axes(axisbg='w',frameon=False)  # bg not visible?
        barlegend = LegendContainer()

        n,bins = degreehisto(self.D['degree'])
        lines = pylab.bar(bins,n,align='center',color='black',width=1.0)
        barlegend.append(lines[0],'total')

        n,bins = degreehisto(self.D['indegree'])
        lines = pylab.bar(bins,n,align='center',color='yellow',width=2./3.)
        barlegend.append(lines[0],'in')

        n,bins = degreehisto(self.D['outdegree'])
        lines = pylab.bar(bins,n,align='center',color='red',width=1./3.)
        barlegend.append(lines[0],'out')

        pylab.legend(*barlegend.args())
        pylab.xlim(xmin=-0.5)
        pylab.xlabel('degree')
        pylab.ylabel('count')
        #pylab.title("Distribution of node degrees (including connections with bulk)")
        pylab.savefig(self.pdf['degree'])

    def _lifetime_histo(self,numfig=1):
        msg(1,"figure %(lifetime)s: life time histogram" % self.pdf)
        pylab.figure(numfig)
        pylab.clf()
        pylab.axes(axisbg='w',frameon=False)  # bg not visible?
        barlegend = LegendContainer()

        n,bins=numpy.histogram(self.D['site_lifetime_avg'],
                               bins=numpy.concatenate((numpy.arange(0,1000,50),
                                                       [1000,5000,10000,20000])))
        xcutoff = 1000
        N_outliers = numpy.sum(n[bins > xcutoff])
        lines = pylab.bar(bins,n,align='edge',color='k',linewidth=0,width=50,
                          label='life time')
        barlegend.append(lines[0],'life time\n(%(N_outliers)d counts > '
                         '%(xcutoff)g)' % locals() )

        pylab.xlim(0,xcutoff)
        pylab.legend(*barlegend.args())
        pylab.xlabel('life time/ps')
        pylab.ylabel('count')
        pylab.savefig(self.pdf['lifetime'])

    def _occupancy_histo(self,numfig=1):
        msg(1,"figure %(occupancy)s: occupancy histograms" % self.pdf)
        pylab.figure(numfig)
        pylab.clf()
        pylab.axes(axisbg='w',frameon=False)  # bg not visible?
        barlegend = LegendContainer()

        n,bins = numpy.histogram(self.D['site_occupancy_kin_avg'],bins=numpy.arange(0.1,1.6,0.1))
        lines = pylab.bar(bins,n,align='edge',color='k',linewidth=0,width=0.1)
        barlegend.append(lines[0],'from kinetics')
        
        n,bins = numpy.histogram(self.D['site_occupancy_rho_avg'],bins=numpy.arange(0.1,1.6,0.1))
        lines = pylab.bar(bins,n,align='edge',color=(0.7,0.7,0.7),alpha=0.7,linewidth=0.1,width=0.1)
        barlegend.append(lines[0],'from density')

        pylab.legend(*barlegend.args())
        pylab.xlabel('occupancy')
        pylab.ylabel('count')
        pylab.savefig(self.pdf['occupancy'])

class LegendContainer:
    """For each bar plot, record first lines instance and the label
    with
    >>> Legend = LegendContainer()
    >>> lines = pylab.bar(...)
    >>> Legend.append(lines[0],'plotlabel')
    Once all legends have been collected, build the legend with
    >>> pylab.legend(*Legend.args())
    """
    def __init__(self):
        self.lines = []
        self.labels = []
    def append(self,line,label):
        self.lines.append(line)
        self.labels.append(label)
    def args(self):
        """Use as pylab.legend(**Legend.args())."""
        return tuple(self.lines), tuple(self.labels)


class HeatmapAnalysis:
    """Combine Hopgraph statistics for a number of simulations into a grid,
    normalize each observable, and color. Clustering is performed if the R
    package is installed in the system. The idea is to quickly compare a number
    simulations based on a combination of observables.
    """
    
    prune_default = ['G_degree_in','G_degree_out',
                     'G_degree_nobulk','G_degree_in_nobulk','G_degree_out_nobulk',
                     'G_degree_min','G_degree_in_min','G_degree_out_min',
                     'G_degree_in_max','G_degree_out_max',
                     'G_degree_in_nobulk_min','G_degree_out_nobulk_min',
                     'G_order_nobulk',
                     'site_N_equivalence_sites','site_N_subsites',
                     ]

    def __init__(self,hoppinggraphs,normalization="maxabs",verbosity=1,prune='default'):
        """Create a 'heatmap' for the Hopgraph statistics from a dictionary of CombinedGraphs.
        
        >>> hm = HeatmapAnalysis(hg,normalize="maxabs")
        
        :Arguments:
        HoppingGraphs    Dictionary of HoppingGraph instances. The key is used to label 
                         the simulation in the heat map and thus should be expressive.          
        normalization    Method to normalize the data across observables. Can be None 
                         (not recommended), 'maxabs', or 'zscore'. See the normalize()
                         method for documentation.
                         NOTE that the normalization strongly influences the clustering 
                         in the heat map.
        verbosity        Chattiness; use at least 1 in order to be notified if you 
                         should install additional packages. Otherwise a less powerful 
                         alternative is chose silently,
        prune            dict with keys that are removed from the heat map; see 
                         prune_default class attribute.

        :Methods:
        plot             plot the heat map
        normalize        normalize using the 'normalize' method
        labels           dictionary of row, column names (and the normalization constants 
                         as strings)
        annotation       'enumerate' dictionaries of labels but not stringified
        """
        # TODO: input a dict of HopGraph or CombinedGraph instances
        #       and a dict with metadata
        # Currently assumes CombinedGraph as input
        set_verbosity(verbosity)

        hg = hoppinggraphs
        self._filename = None   # holds default filename
        
        if prune is 'default':
            prune = self.prune_default
        elif prune is None:
            prune = []
        else:
            try: 'foo' in prune
            except: raise TypeError("prune must support 'in' lookup")

        self.columns = hg.keys()  # labels for sim: 'R13_1_plm', 'R13_1_apo', ...
        self.names = []
        self.heatmap = None
        self.normalization_method = normalization

        # boring data reorganization
        table = {}
        for sim,hopgraph in hg.items():
            for k,v in hopgraph.stats().items():
                try: 
                    table[k].append(v)
                except: 
                    table[k] = [v]
        self.names = [k for k in sorted(table.keys()) if k not in prune]
        self.data  = numpy.array([table[k] for k in self.names])

        # reference values, tagged with row name
        # colum offsets are the indices in self.columns
        self.taggeddata = numpy.rec.fromarrays(self.data,names=self.names)

        self.heatmap = self.data.copy()
        self.normalizations = numpy.zeros(len(self.heatmap))

        self.normalize(method=normalization)  # modifies heatmap and normalizations
        # set up annotation
        self.annotation()

    filename = hop.utilities.filename_function

    __normalization_methods = [None,"maxabs","zscore"]

    def normalize(self,method=None):
        """Normalize the data by row.

        normalize(method=None|'zscore'|'maxabs')

        method can be
        None         Return the unchanged data array.
        'maxabs'     Take the largest absolute value in each row/column and
                     divide each entry in the row/column by it. This results
                     in values between -1 and +1.
        'zscore'     (X-<X>)/sd(X)

        Sets self.heatmap, self.normalizations, self.normalization_method

        normalizations only makes sense for 'maxabs'; in all other cases it
        only contains zeroes.
        """
        if method not in self.__normalization_methods:
            raise NotImplementedError("Normalization '%s' not implemented. Only " % method
                                      +str(self.__normalization_methods)+" are available.")
        self.normalization_method = method
        if method is None:
            self.heatmap = self.data.copy()
            self.normalizations = numpy.zeros(len(self.heatmap))
            return
        # normalize each row separately to max|x|
        self.heatmap = self.data.copy()
        normalizations = []   # holds normalization constant for each row
        for row in self.heatmap:
            if method is "maxabs":
                # crappy but can't figure out the numpy way to set the
                # normalization to the largest absolute value but then
                # store the actual signed value in normalization
                rmin,rmax = row.min(),row.max()
                if abs(rmin) > abs(rmax):
                    rmax = abs(rmin)
                    normalizations.append(rmin)
                else:
                    normalizations.append(rmax)
                if rmax == 0: continue
                row /= rmax
            elif method is "zscore":
                avg = numpy.average(row)
                std = numpy.std(row)
                normalizations.append(avg)
                if numpy.all(std == 0):       # avoid div by zero
                    row *= 0                  # set zscore to 0
                else:
                    row[:] = (row - avg)/std  # z-score
        self.normalizations = numpy.array(normalizations)


    def plot(self,filename=None,format='pdf',**kwargs):
        """Plot the heatmap and save to an image file.

          plot()             # display using windowing system
          plot('hm')         # --> hm.pdf
          plot('hm.png')     # --> hm.png
          plot('hm','png')   # --> hm.png
        
        By default a clustered heat map is constructed using R's heatmap.2
        function. If R cannot be found, an unclustered heat map is
        plotted. **kwargs can be used to customize the output.

        :Arguments:
        filename       name of the image file; may contain extension
                       If empty use the windowing system.
        format         eps,pdf,png... whatever matplotlib understands
        
        **kwargs for R:
        scale          Determines the coloring. Choose between 'none' (the
                       actual values in the heat map (possibly already normalized)),
                       'row' or 'column' (z-score across the dimension)
        N_colors       Number of color levels; default is 32.

        **kwargs for matplotlib:
           The kwargs are applied to the matplotlib.text() method and
           are typically used to set font properties. See the
           pylab/matplotlib documentation.
        """
        if filename:
            format = hop.utilities.fileextension(filename,default=format)
        labels = self.labels()
        try:
            import rpy
            self._heatmap_R(labels,filename=filename,format=format,**kwargs)
        except ImportError:
            msg(0,"rpy package missing: cannot plot clustered heat map, defaulting to "
                "an unclustered heat map")
            self._heatmap_matplotlib(labels,filename=filename,format=format,**kwargs)
        if filename:
            msg(1,"Wrote image to file %s.\n" % self.filename(filename,ext=format))

    def _heatmap_R(self,labels,filename=None,format='pdf',**kwargs):
        """Plot a clustered heat map using R."""
        # Note on coloring: test if data is z-score normalized and if
        #  true, force scale='row' for the heat map so that this is
        #  reflected in the legend (in fact, the z-score is recomputed
        #  over the data for heatmap.2 coloring but it looks identical
        #  and in any case, the clustering is done on the original
        #  data --- which is NOT clear from the heatmap docs)
        from rpy import r, RException
        hm_args = dict(scale='none',
                       margins=(10,10), # space for long labels
                       N_colors=32,
                       )
        hm_args.update(kwargs)
        N_colors = hm_args['N_colors']  # N_colors not a true heatmap argument
        del hm_args['N_colors']
        if filename is not None:
            interactive = False
            def r_format():
                r[format](filename)
            def r_dev_off():
                r.dev_off()
        else:
            interactive = True
            def r_format():
                try:
                    r.X11()
                except RException:
                    r.quartz()
            def r_dev_off():
                msg(1,"Interactive R display. Kill with hop.analysis.kill_R()")
                pass
        if self.normalization_method is 'zscore' and hm_args['scale'] is not 'row':
            hm_args['scale'] = 'row'
            msg(3,"Forcing scale='row' for z-score normalized heat map.\n")
            # (This only has the effect to put the label 'Row Z-score' in the graph)
        try:
            r.library('colorRamps')
            r_color = r.matlab_like(N_colors) # getting somewhat close to matplotlib 'jet'
        except RException:
            msg(1,"For matplotlib-like colors install the R-package 'colorRamps'.")
            r_color = r.topo_colors(N_colors)
        try:
            r.library('gplots')
            r_heatmap = r.heatmap_2
            hm_args.update(dict(key=r.TRUE, symkey=r.FALSE,
                                density_info='histogram', trace='none'))
        except RException:
            msg(1,"For heatmap with a legend install the R-package 'gplots' via \n"
                ">>> import rpy\n"
                ">>> rpy.r.install_packages('gplots',type='source')")
            r_heatmap = r.heatmap
        
        r_format()
        r_heatmap(self.heatmap,
                  labRow=labels['observables'],
                  labCol=labels['columns'],
                  col=r_color,
                  **hm_args)
        r_dev_off()

    def _heatmap_matplotlib(self,labels,filename=None,format='pdf',**kwargs):
        """Plot a un-clustered heat map using matplotlib."""
        import pylab
        pylab.clf()
        pylab.imshow(self.heatmap,interpolation='nearest')
        pylab.axis('off')
        pylab.colorbar(pad=0.075)
        col_labels = [_add_x_ticklabel(n,label,**kwargs) 
                      for n,label in enumerate(labels['columns'])]
        row_labels = [_add_y_ticklabel(n,label,**kwargs) 
                      for n,label in enumerate(labels['observables'])]
        norm_labels = [_add_y_ticklabel(n,label,
                                       offset=self.heatmap.shape[1], # number of colums
                                       horizontalalignment='left',
                                       **kwargs) 
                       for n,label in enumerate(labels['normalizations'])]
        if filename is not None:                
            pylab.savefig(self.filename(filename,ext=format,set_default=True))
        else:
            msg(3,"Only display, no file written because filename == None.")

        
    def annotation(self):
        self.columns_idict = self._make_idict(self.columns)
        self.observables_idict = self._make_idict(self.names)        
        self.normalizations_idict = self._make_idict(self.normalizations)
        return {'columns':self.columns_idict,
                'observables':self.observables_idict,
                'normalizations':self.normalizations_idict}

    def labels(self,precision=2):
        """labels of the columns (simulations) and rows (observables)"""
        return {'columns': map(str,self.columns),
                    'observables': map(str,self.names),
                    'normalizations': ["%g" % round(x,precision) 
                                       for x in self.normalizations]}
    

    def print_annotation(self):
        #pp = pprint.PrettyPrinter() # use global pp
        print "column legend (sims) columns_idict"
        pp.pprint(self.columns_idict)
        print "row legend (observable) observables_idict"
        pp.pprint(self.observables_idict)        

    def _make_idict(self,a):
        return dict(enumerate(a))

def TEST_load_my_sims():
    """Dict of CombinedGraph instances as input for HeatmapAnalysis.
    Hard coded loader for testing.
    """
    ligand = {0:'apo',1:'plm'}  # state --> ligand
    hopgraphs = {}

    # hard coded list for testing, use expressive keys
    # Single hop graphs
    hopgraphs['PME_s_apo'] = hop.graph.HoppingGraph(filename='/Users/oliver/Biop/Projects/WaterNetworks/testcases/1IFC/hopgraph.pickle')
    hopgraphs['PME_s_plm'] = hop.graph.HoppingGraph(filename='/Users/oliver/Biop/Projects/WaterNetworks/testcases/2IFB/hopgraph.pickle')

    # combined hopgraphs. G0 == apo, G1 == holo==plm
    cg = {}
    cg['R13_1'] = hop.graph.CombinedGraph(filename='GSBP/R13/analysis/cg_1IFC_2IFB_1.pickle')
    cg['R15_0'] = hop.graph.CombinedGraph(filename='GSBP/R15/analysis/cg_1IFC_2IFB_0.pickle')
    cg['R15_1'] = hop.graph.CombinedGraph(filename='GSBP/R15/analysis/cg_1IFC_2IFB_1.pickle')
    cg['PME_0'] = hop.graph.CombinedGraph(filename='LAMMPS/100mM/analysis/cg_1IFC_2IFB.pickle')
    cg['PME_TAP']=hop.graph.CombinedGraph(filename='LAMMPS/100mM/analysis/cg_TAP_1IFC_2IFB.pickle')
    cg['PME_CRBP'] = hop.graph.CombinedGraph(filename='../../../CRBPII/comparison/analysis/cg_1OPA_1OPB.pickle')
    
    for sim,combgraph in cg.items():
        for state,hopgraph in enumerate(combgraph.graphs):  # assuming a CombinedGraph
            #print "sim = %(sim)s  state = %(state)d"  % locals()
            sim_id = str(sim)+'_'+ligand[state]
            hopgraphs[sim_id] = hopgraph

    return hopgraphs

def _add_y_ticklabel(y,s,offset=-1,**kwargs):
    kwargs = dict(fontsize=6,
                    verticalalignment='center',horizontalalignment='right')
    textargs.update(kwargs)
    x = offset  # -1 == left of axis, in image coords
    return pylab.text(x,y,s,**textargs)

def _add_x_ticklabel(x,s,offset=-1,**kwargs):
    textargs = dict(fontsize=6,
                    rotation='vertical',
                    verticalalignment='bottom',horizontalalignment='center')
    textargs.update(kwargs)
    y = offset  # -1 == over top of axis, in image coords    
    return pylab.text(x,y,s,**textargs)

def kill_R():
    """Manual last resort to kill the R quartz() window."""
    from rpy import r
    r.dev_off()
    
