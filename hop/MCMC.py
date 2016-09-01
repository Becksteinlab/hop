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
Markov Chain Monte Carlo on hopping graph --- :mod:`hop.MCMC`
=============================================================

The :mod:`hop.MCMC` module uses the information encoded in a hopping
graph to set up a Markov Chain Monte Carlo sampling procedure that
allows one to rapidly generate site occupancy distributions that are
distributed in the same way as the one sampled from MD.

The most convenient entry point is the :func:`hop.MCMC.run`  function ::

  M = MCMC.run(filename='hopgraph.pickle',Ntotal=<int>)

It takes as input a stored hopgraph and immediately runs an MCMC run of *Ntotal*
steps. The output is a :class:`MCMCsampler` object. It contains the 'trajectory' and
useful analysis functions. (Use interactive introspection in ipython to explore
the possibilities of the object.)

Notes on the algorithm:

* some sort of dynamic lattice Monte Carlo with very simple acceptance
  probabilities (0 or 1, if there's no space on the site, and 1 if
  there is)

  ... is 'MCMC' actually the proper description?

* Extension to multiply occupied sites: use the site occupancy
  distributions from siteanalysis, and replace the unconditional move
  by an acceptance probability == s_i(n)

* I am currently using time-forward (out-of) and time-backward (into)
  site moves (the latter inspired by coupling from the past).

Classes and functions
---------------------

"""
from __future__ import absolute_import

import numpy

from . import graph
from .constants import SITELABEL
from . import utilities
from .utilities import msg, set_verbosity


Nequil_default = 10000

class MCMCsampler(utilities.Saveable):
    """Generate an equilibrium distribution of states from a hop graph."""

    _saved_attributes = 'all'

    def __init__(self,h=None,min_hops_observed=1,filename=None):
        """Initialize the Markov Chain Monte Carlo sample with a HoppingGraph.

        M = MCMCsampler(HoppingGraph)

        Build a Markov Chain model from a <HoppingGraph> (with edges deleted
        that have less than <min_hops_observed> hops).
        """
        if h is not None:
            # standard init: setup everything from hopgraph
            import copy
            # 1) compute transition probabilities between all nodes.
            #    represent T_ji (j-->i) as a sparse matrix in dict of dict format
            h.filter(exclude={'Nmin':min_hops_observed})
            self.graph = h.select_graph(use_filtered_graph=True)
            self.site_properties = h.site_properties
            A = self.graph.adj   # adjacency matrix as dict-of-dicts; A[i][j] == j --> i !!
            Nsites = numpy.max(self.graph.nodes())+1
            self.T = numpy.zeros((Nsites,Nsites))   # exit probabilities  T[j,i] = P(i -> j|i)
            self.R = numpy.zeros((Nsites,Nsites))   # entry probabilities R[i,j] = P(j -> i|i?) (?)
            # may look into sparse matrices...

            # fast lookup of neighbours
            self.out_neighbors = dict( [ (i,numpy.sort([e[1] for e in self.graph.out_edges(i)]))
                                         for i in self.graph.nodes()] )
            self.in_neighbors  = dict( [ (i,numpy.sort([e[0] for e in self.graph.in_edges(i)]))
                                         for i in self.graph.nodes()] )

            # probably slow but simple ... actually, fast enough so that I don't care
            # format of a edge in the graph:
            # i : { j1: {'k':k, 'N':N, 'fit':fit_instance}, j2: {...}, ...}
            # i : {
            # {1: {'k':0.100953823261, 'N':32, 'fit':<fitExp2[ 0.95889814  0.10095382  0.00467457]>}},
            # {12: {'k':0.00272870618459, 'N':1, 'fit':<fitExp [0.00272870618459]>}}
            # }

            # probability to exit i to j = k(i->j)/Z with Z = Sum_l k(i-->l)
            # or                           N_ji/Z' with Z' = Sum_l N_li
            # probability to enter i from j:
            #                              N_ij/Y with Y = Sum_l N_il
            for i in self.graph.nodes():
                Nout = numpy.sum( [A[i][j]['N'] for j in self.out_neighbors[i]] )  # total exiting from i
                Nin  = numpy.sum( [A[j][i]['N'] for j in self.in_neighbors[i]] )   # total entering into i
                for j in self.out_neighbors[i]:
                    self.T[j,i] = A[i][j]['N'] * 1.0 / Nout
                for j in self.in_neighbors[i]:
                    self.R[i,j] = A[j][i]['N'] * 1.0 / Nin

            # cumulative prob. distrib. of exit nodes (could sort by decreasing
            # probability to optimize for more probable hits but then must sort
            # nodes as well)
            self.exit_cdf = dict( [(i,numpy.cumsum(self.T[self.out_neighbors[i],i]))
                                   for i in self.graph.nodes() if len(self.out_neighbors[i])>0 ] )
            self.entry_cdf = dict( [(i,numpy.cumsum(self.R[i,self.in_neighbors[i]]))
                                   for i in self.graph.nodes() if len(self.in_neighbors[i])>0 ] )

            self.init_state()

            # obtain maxsite from volume of the site (assume waters are hard spheres)
            # max_occupancy = volume_site/volume_water + 1
            v_water = 4./3.*numpy.pi*hop.constants.r_water**3
            def max_occupancy(sitelabel):
                return int(h.site_properties.volume[sitelabel]/v_water + 1)

            self.maxstate = dict( [(i,max_occupancy(i)) for i in self.graph] )
            self.maxstate[SITELABEL['bulk']] = int(1e10)   # no restrictions on bulk

            # Site probability p_i: probability to find site occupied by at least one particle.
            # p(i) = P(nu_i >= 1) = 1 - P(nu_i = 0)
            # For max nu_i = 1: p(i) = <nu_i>/max nu_i
            # For max nu_i > 1: worry about it later...
            self.siteprobability = dict( [(i,self.site_properties.occupancy_avg[i]/self.maxstate[i])
                                          for i in self.graph] )
            # Set bulk to a random value such as 1/2; it is not being used.
            self.siteprobability[SITELABEL['bulk']] = 0.5
        elif filename is not None:
            # setup from pickled file
            self.load(filename)
        else:
            raise ValueError('Either a hopgraph or a pickled file is needed.')

    def init_state(self,Nbulk=1e4):
        """Initialize state with 1 particle per standard site and Nbulk for the bulk site."""
        # 2) initialize state vector
        #    (many on bulk and 1 on proper sites)
        self.state = dict( [(i,1) for i in self.graph])
        self.state[SITELABEL['bulk']] = int(Nbulk)       # large bulk reservoir
        self.states = numpy.array([])                    # hold 'trajectory'
        self.iteration = 0                               # number of iterations up to/incl states[-1]

    def sample(self,max_iter=10000,record_iterations=True):
        """Run <max_iter> Monte Carlo site moves.

        sample(max_iter=10000)

        Runs a batch of MCMC moves.
        """
        # 3) MCMC sampler
        sites = numpy.array(self.graph.nodes())
        Nsites = len(sites)

        for n in xrange(max_iter):
            site = sites[numpy.random.uniform(SITELABEL['bulk']+1,Nsites)]
            self._MCMCmove(site)
        if record_iterations:
            self.iteration += max_iter

    def run(self,Ntotal=500000,Nskip=1000,verbosity=None):
        """MCMC run multiple cycles of lebgth <Nskip> scans for a total of <Ntotal>.

        run(Ntotal=500000,Nskip=1000)

        Starts from the current configuration in state.
        Creates the collection of configurations states: one state every Nskip steps
        """
        set_verbosity(verbosity)

        Ntotal = int(Ntotal)
        Nskip = int(Nskip)
        Ncycle = int(Ntotal/Nskip)

        self.runparameters = {'Ntotal':Ntotal,'Nskip':Nskip,'Ncycle':Ncycle}
        state_list = []
        msg(1,"Running MCMC:\n\tNtotal = %(Ntotal)d\n\tSaving state every cycle of Nskip steps = %(Nskip)d\n\t"
            "Hence in total Ncycle cycles = %(Ncycle)d\n" % locals())
        for i in xrange(1,Ncycle+1):
            msg(3,"cycle %(i)4d/%(Ncycle)4d\r" % locals())
            self.sample(Nskip)
            state_list.append(self.statevector)
        msg(3,"\n")

        if len(self.states) == 0:
            self.states = numpy.array(state_list)
        else:
            msg(2,"Extending states from %d configurations by %d new ones." \
                    % (len(self.states), len(state_list)))
            self.states = numpy.concatenate((self.states, numpy.array(state_list)))

    def _MCMCmove(self,i):
        if numpy.random.uniform(0,1) < self.siteprobability[i]:
            # try entry move (add one particle to site)
            if self.state[i] >= self.maxstate[i]: # only transfer if there's space on node
                return False
            jnodes = self.in_neighbors[i]       # all entry nodes
            if len(jnodes) == 0:                 # isolated node
                return False
            x = numpy.random.uniform(0,1)        # random number for deciding exit channel
            j = jnodes[ x <= self.entry_cdf[i] ][0]  # j for which p(j-1) < x <= p(j) ('=' necessary to catch x==1)
            if self.state[j] == 0:
                return False       # no particles to transfer from source site
            self.state[i] += 1
            self.state[j] -= 1
            return True
        else:
            # try exit move (remove one particle from site)
            if self.state[i] == 0:
                return False                     # no water to move
            jnodes = self.out_neighbors[i]      # all exit nodes
            if len(jnodes) == 0:                 # isolated node
                return False
            x = numpy.random.uniform(0,1)        # random number for deciding exit channel
            j = jnodes[ x <= self.exit_cdf[i] ][0]  # first jnode for which x > prob ('=' necessary to catch x==1)
            if self.state[j] >= self.maxstate[j]: # only transfer if there's space on node
                return False
            self.state[i] -= 1
            self.state[j] += 1
            return True

    def occupancy(self):
        """Ensemble averaged occupancy (over ALL states incl bulk) and fluctuation."""
        return self.states.mean(axis=0), self.states.std(axis=0)

    def mean(self):
        """Mean for each site (excl bulk)."""
        return self.states[:,self.firstsiteindex:].mean(axis=0)

    def std(self):
        """Standard deviation for each site."""
        return self.states[:,self.firstsiteindex:].std(axis=0)

    def mean_std(self):
        "Returns site labels, mean, and standard deviation for each site (excl bulk)."""
        return self.sites[self.firstsiteindex:], self.mean(), self.std()

    def firstsiteindex():
        doc = """State array index of the first site after bulk."""
        def fget(self):
            return self.node2index[SITELABEL['bulk']]+1
        return locals()
    firstsiteindex = property(**firstsiteindex())

    def statevector():
        doc =  """State as a numpy array; the corresponding nodes are state.keys()"""
        def fget(self):
            return numpy.array(self.state.values())
        def fset(self,x):
            raise NotImplementedError("statevector cannot be assigned; modify state.")
        return locals()
    statevector = property(**statevector())

    def node2index():
        doc = """Translates node label (in graph) to the sequential array index."""
        def fget(self):
            # NOTE: if self.graph is changed manually then this will be wrong
            try:
                return self.__node2index  # use cached dictionary
            except:
                self.__node2index = dict( [(n,idx) for idx,n in enumerate(self.graph.nodes())] )
            return self.__node2index
        return locals()
    node2index = property(**node2index())

    def index2node():
        doc = """Translates sequential array index to node label (in graph)."""
        def fget(self):
            # NOTE: if self.graph is changed manually then this will be wrong
            try:
                return self.__index2node
            except:
                self.__index2node = numpy.array(self.graph.nodes())
            return self.__index2node
        return locals()
    index2node = property(**index2node())
    sites = index2node

    def plot(self,filename=None,plot_skip=None):
        """Plot density plot of the saved configurations in states[]."""
        if len(self.states) == 0:
            raise ValueError("No configurations computed yet. Execute 'run()' first.")

        import pylab

        if plot_skip is None:
            plot_skip = self.runparameters['Nskip']
        plot_step = int(plot_skip/self.runparameters['Nskip']) or 1
        if plot_step == 1:
            plot_skip = self.runparameters['Nskip']

        site_states = self.states[:,self.firstsiteindex:]  # ignore 0 (bulk)
        maxsite = numpy.max([self.maxstate[site] for site in self.graph \
                                 if site <> SITELABEL['bulk']])

        # pylab.clf()
        img = pylab.imshow(site_states[::plot_step].transpose(),
                           interpolation='nearest',cmap=pylab.cm.hot)
        img.axes.set_aspect('auto')
        pylab.colorbar(ticks=range(maxsite+1),shrink=0.75)

        pylab.title('Number of water molecules per site')
        pylab.xlabel('MCMC iteration/%d' % plot_skip)
        pylab.ylabel('site (sequentially enumerated)')

        pylab.draw_if_interactive()
        if filename:
            pylab.savefig(filename)

    def plot_occupancy(self,legend=True,**kwargs):
        """Plot site label vs <N> +/- std(N).

        legend      True: add legend, False: return (line,description)
        **kwargs    additional arguments for errbar plot such as color='k', fmt='o'

        """
        import pylab
        plotargs = {'fmt':'ko', 'color':'black'}
        plotargs.update(kwargs)
        sites = self.index2node
        o_mean, o_std = self.occupancy()
        firstsiteindx = self.node2index[SITELABEL['bulk']]+1
        # upper graph: mean + stdev bar
        ##pylab.subplot(211,frameon=False)
        pylab.subplot(111,frameon=False)
        pylab.title('occupancy' )
        pylab.ylabel(r'$<n_i>$')
        lines = pylab.errorbar(x=sites[firstsiteindx:],
                               y=o_mean[firstsiteindx:],
                               yerr=o_std[firstsiteindx:],
                               **plotargs)
        labelstring = ''   # no legend
        #labelstring = r'$N=%g$' % self.runparameters['Ntotal']
        #if legend:
        #    pylab.legend((lines[0],),(labelstring,),
        #                 loc='upper right',numpoints=1,shadow=True)
        pylab.xlabel('site label')

        # lower graph: stdev
        ##pylab.subplot(212,frameon=False)
        ##pylab.xlabel('site label')
        ##pylab.ylabel(r'$\sigma(\nu_i)$')
        ##pylab.plot(sites[firstsiteindx:],
        ##           o_std[firstsiteindx:],'o',
        ##           color=plotargs['color'])

        if not legend:
            return (lines[0],labelstring)
        return None

    def occupancy_mean_correl(self):
        """Calculates the correlation coefficient between simulation and MCMC occupancies."""
        corr = numpy.corrcoef(*self._occupancies_mean())
        return corr[0,1]

    def occupancy_std_correl(self):
        """Calculates the correlation coefficient between simulation and MCMC occupancy fluctuations."""
        corr = numpy.corrcoef(*self._occupancies_std())
        return corr[0,1]

    def plot_correl(self,legend=True,**kwargs):
        """Plot the occupancy from the MD simulation vs the MCMC one."""
        import pylab

        occ_MD,occ_MCMC = self._occupancies_mean()
        corr = self.occupancy_mean_correl()
        labelstring = r'$C=%g$' % (corr)
        lines = pylab.plot(occ_MD, occ_MCMC, 'o', label=labelstring, **kwargs)

        # straight line fit
        import graph
        f = graph.fitlin(occ_MD,occ_MCMC)
        x = numpy.array([0,1])
        pylab.plot(x,f.fit(x),'--',**kwargs)

        pylab.xlabel('MD occupancy')
        pylab.ylabel('MCMC occupancy')
        if legend:
            pylab.legend(loc='best',numpoints=1,
                         prop=pylab.matplotlib.font_manager.FontProperties(size='smaller'))
        if not legend:
            return (lines[0],labelstring)
        return None

    def _occupancies_mean(self):
        """Returns the mean occupancy from MD and MCMC."""
        occ_MD = self.site_properties.occupancy_avg[self.index2node]
        occ_MCMC = self.occupancy()[0]
        firstsiteindx = self.node2index[SITELABEL['bulk']]+1
        return occ_MD[firstsiteindx:], occ_MCMC[firstsiteindx:]

    def _occupancies_std(self):
        """Returns the fluctuation (stdev) in the occupancy from MD and MCMC."""

        raise NotImplementedError('MD Occupancy must be read from SiteAnalysisArray.')

        occ_MD = self.site_properties.occupancy_std[self.index2node]
        occ_MCMC = self.occupancy()[1]
        firstsiteindx = self.node2index[SITELABEL['bulk']]+1
        return occ_MD[firstsiteindx:], occ_MCMC[firstsiteindx:]

    def autocorrelation(self,start=None,stop=None,step=None,**kwargs):
        """Calculates the auto correlation function for all site trajectories."""
        from .utilities import autocorrelation_fft as ACF
        return numpy.array([ACF(traj,**kwargs)
                              for traj in self.states.T[:,start:stop:step]])

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
        tmp = numpy.array([avACF(traj,**kwargs) for traj in self.states.T[:,::step]])
        return tmp[:,0], tmp[:,1]


def run(filename='hopgraph.pickle',Ntotal=500000,Nskip=1000,Nequil=Nequil_default):
    """Perform Markov Chain Monte Carlo on a model derived from the hopping graph."""

    h = graph.HoppingGraph(filename=filename)
    M = MCMCsampler(h)

    msg(0,"MCMCsampler() for %s\n" % filename)

    if Nequil > 0:
        # burn in/flooding
        msg(1,"Running %(Nequil)d steps for equilibration...\n" % locals())
        M.sample(Nequil,record_iterations=False)

    # MCMC run
    M.run(Ntotal=Ntotal,Nskip=Nskip,verbosity=3)

    return M

class Pscan(object):
    """Run a MCMC sampler for a number of parameter values."""

    def __init__(self,parameter,pvalues=None,filename='hopgraph.pickle',
                 Ntotal=1e6,**kwargs):
        """Sample on hopping graph for different values of <parameter> = p.

        P = Pscan(parameter=<string>,pvalues=<sequence>,filename=<filename>,**kwargs)

        <parameter> must be a keyword argument to hop.MCMC.run();
        the parameter overrides any default values that may have been
        set. For instance, <parameter> can be 'Ntotal' or 'filename'.

        kwargs: all other kwargs are directly passed on to MCMC.run().
        """

        self.parameter = parameter
        self.pvalues = list(pvalues)
        self.Ntotal = Ntotal
        M = {}
        for p in pvalues:
            print "parameter %s = %g" % (parameter,p)
            kwargs[parameter]=p
            kwargs.setdefault('Ntotal',Ntotal)
            kwargs.setdefault('filename',filename)
            M[p] = run(**kwargs)

        self.MCMCsamplers = M
        del M

    def occupancy_mean_correl(self):
        """Returns X=pvalues, Y=occupancy_mean_correlations."""
        self.pvalues.sort()
        return (numpy.array(self.pvalues),
                numpy.array([self.MCMCsamplers[v].occupancy_mean_correl()
                             for v in self.pvalues]))

    def plot_occupancy_mean_correl(self,linestyle='ko-',**kwargs):
        import pylab
        x,y = self.occupancy_mean_correl()
        lines = pylab.semilogx(x,y,linestyle,**kwargs)
        pylab.xlabel('total MCMC steps')
        pylab.ylabel('correlation coefficient with MD')

    def plot_states(self,maxcolumns=2):
        """Plot all state 'trajectories' as a tiled plot."""
        import pylab
        pvalues = sorted(self.MCMCsamplers.keys())
        Nmax = len(pvalues)
        maxrows = int(Nmax*1.0 / maxcolumns + 0.5)
        nplot = 0
        pylab.clf()
        for irow in xrange(maxrows):
            row = irow + 1
            for icolumn in xrange(maxcolumns):
                column = icolumn + 1
                nplot += 1
                iplot = nplot - 1
                p = pvalues[iplot]
                pylab.subplot(maxrows,maxcolumns,nplot, frameon=False)
                self.MCMCsamplers[p].plot()
                ax = pylab.gca()
                if ax.is_first_col():
                    pylab.ylabel('site #')
                else:
                    pylab.ylabel('')
                if ax.is_last_row():
                    pass # keep the original label
                else:
                    pylab.xlabel('')
                pylab.title(r'%s=%g' % (self.parameter,p))


    def plot_occupancy(self,**kwargs):
        """Plot <n_i> (site occupancy from MCMC) for all parameter values.

        (See _plotter())"""
        self._plotter('plot_occupancy',fmt='o',**kwargs)

    def plot_occupancy_mean_correl(self,**kwargs):
        """Plot MD occupancy vs MCMC occupancy.

        plot_correl(colorscale='linear'|'log')

        (See _plotter())"""
        import pylab
        self._plotter('plot_correl',**kwargs)
        pylab.plot([0,1],[0,1],'k-.')

    def _plotter(self,func,colorscale='linear',**fargs):
        """Plot func for all MCMCsamplers in the collection.

        _plotter('<functionname>',colorscale='linear'|'log',**fargs)

        '<functionname>'      string, names a plotting method such as 'plot_correl'
        colorscale            scale color based on the parameter value ('linear' or 'log';
                              note that 'log' can only be used with numerical data > 0)
        **fargs               additional arguments that the function may need
        """
        import pylab
        # currently, pvalues must be numerical; if not then one should use the enumeration
        # instead of the pvalues list
        if numpy.all(numpy.isreal(self.pvalues)):
            p_dict = dict( zip(self.pvalues, map(float,self.pvalues)) )
        else:
            p_dict = dict( [(p,n) for n,p in enumerate(self.pvalues)] )

        _scalings = {'linear': lambda x:x,
                     'log': lambda x:numpy.log(x)}
        try:
            scale = _scalings[colorscale]
        except KeyError:
            raise ValueError("<colorscale> must be one of "+str(_scalings.keys()))
        pnormalize = pylab.normalize(vmin=numpy.min(scale(p_dict.values())),
                                     vmax=numpy.max(scale(p_dict.values())))

        pylab.clf()
        ax = pylab.axes(frameon=False)

        lines = []
        descriptions = []
        for p in sorted(self.MCMCsamplers):
            p_num = p_dict[p]
            plot_func = self.MCMCsamplers[p].__getattribute__(func)
            (line,correl_legend) = plot_func(color=pylab.cm.jet(pnormalize(scale(p_num))),
                                         legend=False, **fargs)
            descr = r'$%s\ %s$' % (str(p), correl_legend.replace(r'$',''))   # minimalistic legend
            lines.append(line)
            descriptions.append(descr)

        pylab.legend(lines,descriptions,numpoints=1,shadow=True,
                     prop=pylab.matplotlib.font_manager.FontProperties(size='smaller'))

    def save(self,filename='pscan.pickle'):
        """Save pscan object to pickle file.

        save(pscan.pickle)

        Load with

            import cPickle
            myPscan = cPickle.load(open('pscan.pickle'))
        """
        import cPickle
        fh = open(filename,'wb')
        try:
            cPickle.dump(self,fh,cPickle.HIGHEST_PROTOCOL)
        finally:
            fh.close()
        msg(3,"Wrote Pscan object to '%s'.\n" % filename)

class MultiPscan(utilities.Saveable):
    """Run Pscan(**pscanargs) <repeat> times and collect all Pscan objects in list.
    """

    _saved_attributes = 'all'

    def __init__(self,repeat=10,**pscanargs):
        """pscans = MultiPscan(repeat=10, parameter='Ntotal', pvalues=[1e4,2.5e4,....], ...)
        See Pscan() for description of pscanargs.
        """
        self.collection = [Pscan(**pscanargs) for i in xrange(repeat)]

    def plot(self,**kwargs):
        multi_plot(self.collection,plottype='whisker',funcname='occupancy_mean_correl',**kwargs)

def multi_plot(plist,plottype='whisker',Nequil=Nequil_default,funcname='occupancy_mean_correl',**kwargs):
    """Display a collection of functions.

    multi_plot(plist,plottype='whisker',Nequil=10000,funcname='occupancy_mean_correl',**kwargs)

    The function is obtained from a method call on the objects in plist. The
    assumption is that these are functions of Ntotal (if not, set Nequil=0; Nequil is
    added to x). Each object is a different realization, e.g. multiple MCMC runs.

    plottype       'whisker' (whisker plot), 'standard' (average and standard deviations)
    Nequil         correction, added to x
    funcname       string; a method of the objects in plist that does EXACTLY the following:
                   x,y = obj.funcname()
                   where x and y are numpy arrays of equal length
    **kwargs       color, boxcolor, mediancolor, capsize
    """
    import pylab

    plottypes = ('whisker','standard')

    # sanity checks
    if plottype not in plottypes:
        raise ValueError("Only plottypes from %(plottypes)r, not %(plottype)s." % locals())

    try:
        plist[0].__getattribute__(funcname)()   # How do I handle old-style classes? Check __getattr__?
    except:
        raise ValueError("funcname='%(funcname)r' has the wrong signature or is not a method "
                         "of the objects in plist." % locals())

    def xy_from(obj):
        return obj.__getattribute__(funcname)()
    def x_from(obj):
        return xy_from(obj)[0]
    def y_from(obj):
        return xy_from(obj)[1]

    kwargs.setdefault('color','b')           # main color
    kwargs.setdefault('capsize',0)

    x = x_from(plist[0]) + Nequil
    all_y = numpy.array([y_from(p) for p in plist])
    ymean = all_y.mean(axis=0)
    yerr = all_y.std(axis=0)

    ax = pylab.axes(frameon=False)
    if plottype == 'standard':
        lines = pylab.errorbar(x,ymean,yerr=yerr,fmt='o',linestyle='-',**kwargs)

        for p in plist:
            px,py = xy_from(p)
            px += Nequil
            components = pylab.semilogx(px,py,'o',color=kwargs['color'],ms=3)
    elif plottype == 'whisker':
        # widths that look the same in a log plot
        def logwidth(x,delta):
            """Return linear widths at x that, logarithmically scaled, appear to be delta wide."""
            q = numpy.exp(delta)
            return x * (1 - q)/(1 + q)

        kwargs.setdefault('mediancolor','k')     # for median in whiskerplot
        kwargs.setdefault('boxcolor','darkgray')

        pylab.semilogx(x,ymean,'o-',color=kwargs['color'],ms=3)
        components = pylab.boxplot(all_y,positions=x,
                                   widths=logwidth(x,kwargs.setdefault('widths',0.15)))
        ax = pylab.gca()
        ax.set_xscale('log')
        ax.set_yscale('linear')
        # must modify components for customisation
        for l in components['boxes']:
            l.set_color(kwargs['boxcolor'])
        for l in components['whiskers']:
            l.set_color(kwargs['color'])
            l.set_linestyle('-')
        for l in components['caps']:
            l.set_color(kwargs['color'])
        for l in components['medians']:
            l.set_color(kwargs['mediancolor'])
            l.set_linewidth(3.0)
        for l in components['fliers']:
            l.set_color(kwargs['color'])

        pylab.draw_if_interactive()

    pylab.xlabel('total MCMC steps')
    pylab.ylabel('correlation coefficient with MD')
    pylab.title('convergence of occupancy correlation with MD')

    #return x,ymean,yerr
