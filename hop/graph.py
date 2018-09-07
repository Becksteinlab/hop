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
Generating and analyzing a hopping graph --- :mod:`hop.graph`
=============================================================

Interprete the high density sites as graph ('transport graph'), with
the sites as vertices and transitions (sampled by the simulation) as
edges. The graph is directed.

Each edge (transition) is decorated with the dominant transition rate,
the number of events seen, and an instance of fit_func, which
represents the fitted function to the survival times.

Each vertex (site) is decorated with the average residency time (and stdev, N).

Typical use of the module::

  TN = TransportNetwork(hoppingTrajectory,density)
  hopgraph = TN.HoppingGraph()
  hopgraph.save('hopgraph')

The basic object is the :class:`hop.graph.HoppingGraph`; see its
documentation for further analysis methods.

Classes and functions
---------------------
"""
from __future__ import absolute_import

import sys
import cPickle
import os.path
import warnings
from itertools import izip

import networkx as NX
import numpy
import scipy.optimize
import Bio.PDB
from MDAnalysis.lib.log import ProgressMeter

import hop
from . import constants
from . import sitemap
from . import trajectory
from .constants import SITELABEL
from .exceptions import MissingDataError, MissingDataWarning
from . import utilities
from .utilities import msg,set_verbosity, iterable, asiterable, CustomProgressMeter




import logging
logger = logging.getLogger("MDAnalysis.analysis.hop.graph")

class TransportNetwork(object):
    """A framework for computing graphs from hopping trajectories.

    The unit of time is ps.

    The :class:`TransportNetwork` is an intermediate data structure
    that is mainly used in order to build a :class:`HoppingGraph` with
    the :meth:`TransportNetwork.HoppingGraph` method.
    """

    def __init__(self, traj, density=None, sitelabels=None):
        """Setup a transport graph from a hopping trajectory instance.

        ::
            hops = hop.trajectory.HoppingTrajectory(hopdcd='whop.dcd',hoppsf='whop.psf')
            tn = TransportNetwork(hops)
        """
        self._cache = dict()
        if not isinstance(traj, trajectory.HoppingTrajectory):
            errmsg = "'traj' must be a <hop.trajectory.HoppingTrajectory> instance."
            logger.fatal(errmsg)
            raise TypeError(errmsg)
        self.traj = traj
        self.n_atoms = self.traj.hoptraj.n_atoms
        self.dt = self.traj.hoptraj.dt  # ps
        self.totaltime = self.traj.totaltime
        self.graph = NX.DiGraph(name='Transport network between sites')

        if not density is None and not isinstance(density, sitemap.Density):
            raise TypeError('density must be a <hop.sitemap.Density> instance.')
        self.density = density  # only needed to pass down site_properties to hopgraph

        if sitelabels:
            self._cache['sitelabels'] = sitelabels

    def graph_alltransitions(self):
        """Constructs the graph that contains all transitions in the trajectory.

        Populates TransportGraph.graph with a graph that contains all
        sites and one edge for each transition that was observed in
        the trajectory. Useful for an initial appraisal of the
        complexity of the problem.

        .. Warning:: Erases any previous contents of graph.
        """
        self.graph = NX.DiGraph(name="Transitions between all sites, including interstitial and outliers")

        pm = ProgressMeter(self.traj.n_frames, interval=100,
                           format="Analyzing transitions: frame %(step)5d/%(numsteps)d  [%(percentage)5.1f%%]\r")


        self.traj.hoptraj.rewind()
        hops = self.traj.iterator()

        # first frame is only needed for slast
        ts = hops.next()
        s = ts._x
        slast = s.astype(int).copy()   # copy of the last frame
        for ts in hops:                # continue with second frame
            pm.echo(ts.frame)
            s = ts._x.astype(int)
            self.graph.add_edges_from(izip(slast, s))
            slast = s.copy()           # copy of the last frame
        self._cache['sitelabels'] = sorted(self.graph.nodes())

    def HoppingGraph(self,verbosity=3):
        """Compute the HoppingGraph from the data and return it."""
        self.compute_residency_times(verbosity)
        return self.hopgraph

    def compute_residency_times(self,verbosity=3):
        self._analyze_hops(verbosity=verbosity)
        logger.info("Building hopping graph with rate constants k1, k2 or k.")
        logger.info("(Ignore 'Warning: ... maxfev = 800': then a single k is chosen automatically.)")
        try:
            site_properties = self.density.site_properties
        except AttributeError:
            site_properties = None
        # should  add more trjdata attr: hoppsf,hopdcd,density(file)
        self.hopgraph = HoppingGraph(self.graph,self.graph_properties,
                                     trjdata={'dt':self.dt,'totaltime': self.totaltime,
                                              'time_unit':'ps'},
                                     site_properties=site_properties)

    def _analyze_hops(self,verbosity=0):
        """Core function to compute the waiting times tau_ji for hops i->j

        :Results:

        self.graph
            directed graph of all transitions
        self.graph_properties
            Each edge (i,j) has a dictionary associated, containing
            the waiting times tau for hops from i --> j, the barrier
            crossing times tbarrier (the time the particle was in the
            interstitial when hopping), the atom index numbers iatom,
            and the frames when site i was exited for the hop.

            iatom can be translated to the original atom serial
            numbers with the hopping trajectory psf.

        (Note: runs a few percent faster with verbosity=0)
        """

        #set_verbosity(verbosity) # not use with ProgressMeter
        try:
            if len(self._cache['sitelabels']) < 3:   # 3 = outlier + 2
                errmsg = 'TransportNetwork: Need at least 2 sites to analyze transitions: N=%d' % \
                    len(self._cache['sitelabels'])
                logger.error(errmsg)
                raise ValueError(errmsg)
        except KeyError:
            pass

        self.graph = NX.DiGraph()
        self.graph.name = "Transitions between sites"
        self.graph_properties = dict()   # contains dicts for nodes and edges

        self.traj.hoptraj.rewind()
        hops = self.traj.iterator()

        s = self.traj.ts._x.astype(int)    # much faster than  numpy.array(snow,dtype=numpy.Int)
        slast = s.copy()                   # unavoidable copy of the last frame
        hops.next()                        # move to second frame

        # state[i] = dict( s0, t0, t1) + slast, s
        #  s0      assigned site        t0      time/frame for entering s0
        #  last    site at t-1          t1      exit from s0
        #  now     site at t           (t2      entering new site != s0)
        # last == slast
        # now  == s
        # (A flat natoms x 5 int array may be faster but for the time being I need
        # clean-ish code...)
        # (could use a numpy.rec --- would do exactly what I want)
        state = numpy.empty(self.n_atoms,dtype=dict)
        state[:] = [{'s0':s[iatom], 't0': self.traj.ts.frame, 't1': None}
                    for iatom in xrange(self.n_atoms)]
        pm = ProgressMeter(self.traj.n_frames, interval=100,
                           format="Analyzing hops: frame %(step)5d/%(numsteps)d  [%(percentage)5.1f%%]\r")

        for ts in hops:
            pm.echo(ts.frame)
            # start from frame 2 (as we already used hops.next() above)
            s = self.traj.ts._x.astype(int)

            # Bugs of the main loop:
            #  * still counts (-1,1) transitions (outlier -> bulk); currently they are
            #    explicitly filtered when building the CombinedGraph
            #  * Dan Willenbring reports that also (-1, N) show up.
            #  * code is ugly and not optimized: should be possible to do this on all atoms at once
            for iatom in xrange(self.n_atoms):
                if s[iatom] == SITELABEL['outlier']:
                    # outliers: count them as whence they came (typically, interstitial or bulk)
                    s[iatom] = slast[iatom]

                if s[iatom] != state[iatom]['s0']:
                    if slast[iatom] == state[iatom]['s0']:
                        # particle exited site s0, so remember the exit time for the future:
                        state[iatom]['t1'] = ts.frame
                    if s[iatom] != SITELABEL['interstitial']:
                        # particle enters a new site s
                        if state[iatom]['s0'] != SITELABEL['interstitial']:
                            # s0 != interstitial ignores trajectories starting off-site,
                            # for all others we collect data now, using the last known exit time:
                            edge = (state[iatom]['s0'], s[iatom])
                            tau = state[iatom]['t1'] - state[iatom]['t0']
                            tbarrier = ts.frame - state[iatom]['t1']
                            try:
                                self.graph_properties[edge]['tau'].append(tau)
                                self.graph_properties[edge]['tbarrier'].append(tbarrier)
                                self.graph_properties[edge]['frame'].append(ts.frame)
                                self.graph_properties[edge]['iatom'].append(iatom)
                            except KeyError:
                                self.graph.add_edge(*edge)
                                self.graph_properties[edge] = {'tau':[tau],'tbarrier':[tbarrier],
                                                               'frame':[ts.frame],'iatom':[iatom]}
                        #else:
                            # particle started from interstitial and it's the first time
                            # on a site so we only re-start its clock ('reset state')
                            # pass
                        # reset state
                        state[iatom] = {'s0':s[iatom], 't0': ts.frame, 't1': None}
            slast = s.copy()   # copy of the last frame

        # mop up last frame: we can't compute rates but at least use
        # the residency time for anything that still sits on a site or
        # sat on a site and is now back in the interstitial
        for iatom in xrange(self.n_atoms):
            if state[iatom]['s0'] == SITELABEL['interstitial']:
                continue
            node = state[iatom]['s0']          # last site visited
            frame = self.traj.ts.frame
            tau = frame - state[iatom]['t0']   # time since ENTERING the last site
            tbarrier = None
            try:
                self.graph_properties[node]['tau'].append(tau)
                self.graph_properties[node]['tbarrier'].append(tbarrier)
                self.graph_properties[node]['frame'].append(frame)
                self.graph_properties[node]['iatom'].append(iatom)
            except KeyError:
                self.graph.add_node(node)
                self.graph_properties[node] = {'tau':[tau],'tbarrier':[None],
                                               'frame':[frame],'iatom':[iatom]}

        # convert times from frames to ps and use numpy arrays to organize data
        # (less clear but easier to slice and dice)
        for n,d in self.graph_properties.items():
            d['frame'] = numpy.array(d['frame'])
            d['iatom'] = numpy.array(d['iatom'])
            d['tau'] = self.dt * numpy.array(d['tau'])
            try:
                d['tbarrier'] = self.dt * numpy.array(d['tbarrier'])
            except TypeError:
                d['tbarrier'] = numpy.array(d['tbarrier'])
            # axis 0 in array is the event number: [tau,tbarrier,iatom,frame]
        if 'sitelabels' not in self._cache:          # update cache
            self._cache['sitelabels'] = sorted(self.graph.nodes())

    def compute_site_times(self,verbosity=3):
        """Compute the 'residency' time of each water molecule on each site.

        The 'site time' of a site i is computed as::

          theta[i] = 1/T_sim ( Sum_j tau[j,i] + Sum tau[i] )

        tau[j,i] is the waiting time for hops from i to j. tau[i] is
        the waiting time for molecules that are not observed to leave
        site i.


        The function updates self.theta[site] for each site with an
        array of residency times (in ps).

        It uses the residency times and thus requires
        compute_residency_times() was run previously.

        :TODO:
          Maybe use the barrier time as well (or a portion thereof,
          perhaps proportional to the barrier height (related to the
          kji) --- rate theory??)
        """
        try:
            props = self.hopgraph.properties
        except AttributeError:
            raise MissingDataError("No hopgraph found. Run 'compute_residency_times()' first.")
        self.theta = dict()        # dict of arrays
        for edge in props:
            if type(edge) is tuple:
                i,j = edge       # i --> j
            else:
                i = edge         # molecule did not leave site
            # could look at 'tbarrier' as well
            try:
                self.theta[i] = numpy.append(self.theta[i], props[edge]['tau'])
            except KeyError:
                self.theta[i] = props[edge]['tau'].copy()       # make a new array

    def compute_site_occupancy(self):
        """Computes occupancies from the residency times theta and updates self.occupancy.

        occupancy::
                      N_i
           o[i] = 1/T Sum theta[i,k]
                      k=1

        where T is the total trajectory time and the sum runs over
        all residency times that were recorded for the site i.

        :Attributes:

        self.occupancy
                         numpy array with occupancies, site label == index
        self.occupancy_error
                         numpy array with error estimates for
                         occupancies (Delta = Delta(theta)/T; this is
                         a biased estimate because Delta(theta) is
                         calculated with N instead of N-1)
        """
        try:
            if len(self.theta) == 0:
                raise ValueError
        except (AttributeError,ValueError):
            # need to compute thetas first
            self.compute_site_times(verbosity=0)

        # array for accumulated site times
        # assumptions: labels are integers from 0 to smax
        smax = numpy.max(self._cache['sitelabels']) # exists after compute_site_times()
        self.occupancy = numpy.zeros(smax+1)        # include 0 so that label matches index
        self.occupancy_error = numpy.zeros(smax+1)  # error estimate (biased, N, not N-1)
        for label,thetas in self.theta.items():
            self.occupancy[label] = numpy.sum(thetas) / self.traj.totaltime
            self.occupancy_error[label] = numpy.std(thetas) / self.traj.totaltime # biased

    def plot_residency_times(self,filename,bins=None,exclude_outliers=True):
        """Plot histograms of all sites

        plot_residency_times('sitetime.eps')

        pylab always writes the figure to the named file. If pylab is
        already running, display the graph with pylab.show().

        The histograms are normalized and the time values are the left
        edges of the bins.

        If bins=None then the number of bins is determined heuristically.
        """
        import pylab

        if exclude_outliers:
            firstSitePos = SITELABEL['bulk']
        else:
            firstSitePos = SITELABEL['interstitial']
        if bins:
            def nbins(data,bins=bins):
                return bins
        else:
            def nbins(data):
                from math import floor,ceil

                nbins = ceil(len(data)/10.0)  # uniform dist: 10 counts per bins
                minresolution = self.dt
                trange = abs(numpy.max(data)-numpy.min(data))
                resolution = trange/nbins

                if nbins < 5: nbins = 5
                if resolution < minresolution:
                    nbins = floor(trange/minresolution)
                return nbins

        pylab.clf()
        for site in self._cache['sitelabels'][firstSitePos:]:
            h,e = numpy.histogram(self.theta[site],
                                  bins=nbins(self.theta[site]),normed=True)
            pylab.plot(e,h, label='site '+str(site), linewidth=2.0)
        pylab.xlabel('residency time  theta/ps')
        pylab.legend()
        pylab.savefig(filename)

    def plot_site_occupancy(self,filename,bins=10,
                            exclude_sites=[SITELABEL['interstitial'],SITELABEL['bulk']]):
        """Plot site occupancy (from compute_site_occupancy().

        plot_site_occupancy(filename,exclude_sites=[0,1])

        filename         name of output file
        bins             bins for histogram (see numpy.histogram)
        exclude_site     list of site labels which are NOT plotted. Typically,
                         exclude interstitial and bulk.
        """
        import pylab

        try:
            labels = range(len(self.occupancy))
        except AttributeError:
            self.compute_site_occupancy()
            labels = range(len(self.occupancy))
        try:
            map(labels.remove, exclude_sites)
        except ValueError:
            pass      # BUG: fails on first non-existent site but should be ok
        labels = numpy.array(labels)
        occ = self.occupancy[labels]
        occ_err = self.occupancy_error[labels]

        h,e = numpy.histogram(occ,bins=bins,normed=True)

        edges = numpy.append(e,[2*e[-1] - e[-2]]) # guess last upper edge
        midpoints = 0.5*(edges[1:]+edges[:-1])

        pylab.clf()
        pylab.subplot(1,2,1, frameon=False)
        pylab.errorbar(labels,occ,occ_err,fmt='k.',capsize=0,linewidth=2.0)
        pylab.xlabel('site label')
        pylab.ylabel('site occupancy')
        pylab.subplot(1,2,2, frameon=False)
        pylab.plot(midpoints,h,'k-',linewidth=2.0)
        pylab.xlabel('site occupancy')
        pylab.legend(('normalized\ndistribution',),'best')

        pylab.savefig(filename)


    def export(self,filename=None,format='png',exclude_outliers=False):
        """Export the graph as dot file and as an image.

        export(filename)

        See https://networkx.lanl.gov/reference/pygraphviz/pygraphviz.agraph.AGraph-class.html#draw for possible output formats.

        See http://graphviz.org/doc/info/attrs.html for attributes.

        Note: On Mac OS X 10.3.9+fink the pygraphviz rendering is
        buggy and does not include node labels. Simply use the
        exported .dot file and use Mac OS X graphviz from
        http://www.pixelglow.com/graphviz/
        """

        import pygraphviz

        gfxname = self.filename(filename,format)
        dotname = self.filename(filename,'dot')       # load into graphviz

        G = pygraphviz.AGraph(directed=True)
        G.add_nodes_from(self.graph.nodes())
        G.add_edges_from(self.graph.edges())

        G.graph_attr['label'] = 'Water transport network in 1IFC'
        G.graph_attr['splines'] = 'true'
        G.graph_attr['fontcolor'] = 'white'
        G.node_attr['shape'] = 'circle'
        G.node_attr['color'] = 'black'
        G.edge_attr['color'] = 'black'
        G.edge_attr['dir'] = 'forward'

        G.layout(prog='circo')
        G.write(dotname)
        G.draw(path=gfxname,format=format)  # buggy on Mac OS X


class HoppingGraph(object):
    """A directed graph that describes the average movement of
    molecules between different well-defined sites by treating the
    sites as nodes and transitions as edges.

    :Attributes:
    graph
        graph with edges; edges contain rates, fit functions, etc
    properties
        raw data for edges
    trjdata
        metadata of the original trajectory
    site_properties
        density-derived node properties, imported from hop.sitemap.Density
    theta
        dict of nodes with residence times (see compute_site_times())
    occupancy_avg
        average occupancy with standard deviation (see compute_site_occupancy())
    occupancy_std
        (numpy array)

    :Methods:
    compute_site_occupancy()
        Computes occupancies from the residency times theta and
        updates self.occupancy_avg and self.occupancy_std.
    compute_site_times()
        Computes residency time theta.
    save()
        save graph as a pickled file
    load()
        reinstantiate graph from saved file; typically just use the
        constructor with the filename argument
    filter()
        make a filtered graph for further analysis and visualization;
        most plot/export functions require a filtered graph
    plot_fits()
        plot fits of the survival time against the data
    tabulate_k()
        table of rate constants
    export()
        export graph as a dot file that can be used with graphviz
    export3D()
        export graph as a psf/pdb file combination for visualization in VMD

    Properties for nodes are always stored as numpy arrays so that one
    can directly index with the node label (==site label), which is an
    integer from 0 to the number of nodes. Note that 0 is the
    interstitial (and only contains bogus data or None), and 1 is the
    bulk. The bulk site is often excluded from analysis because it is
    different in nature from the 'real' sites defined as high density
    regions.
    """

    def __init__(self,graph=None,properties=None,filename=None,trjdata=None,
                 site_properties=None):
        """Directed graph with edges containing the rate k_ji,number of observations and S(t) fit.

          h = HoppingGraph(graph,properties)
          h = HoppingGraph(filename='HoppingGraph.pickle')

        :Arguments:
        graph
           networkx graph with nodes (i) and edges (i,j)
        properties
           dictionary of edges: For each edge e, properties contains a
           dictionary, which contains under the key 'tau' a list of
           observed waiting times tau_ji.  nodes are also listed if
           they do not participate in a transition
        trjdata
           dictionary describing properties of the trajectory
           such as time step 'dt' or name of 'dcd' and 'psf'.

           Attributes that are in use:
              dt
                 time between saved snapshots in ps
              hoppsf
                 hopping trajectory psf file name
              hopdcd
                 hopping trajectory dcd file name
              density
                 pickle file of the density with the sites
              totaltime
                 length of trajectory in ps[*]_

           Not used:
              time_unit
                 'ps'

           .. _[*]: required for :meth:`compute_occupancy`

        site_properties
           list of site properties:
           :attr:`hop.sitemap.Density.site_properties` (add if you want graphs
           with mapped labels) (**Really required for most things...!**)

        When the graph is built from edges and properties then the
        rate constants are calculated. For graphs with many hopping
        events this can take a long time (hours...).

        The decorated and directed graph is accessible as :attr:`HoppingGraph.graph`

        :BUGS:
        * *trjdata* is  required for full functionality but it is currently the user's
          responsibility to fill it appropriately (although :meth:`TransportNetwork.compute_residency_times`
          already adds some data)
        * *site_properties* **are** required and must be added with the constructor
        """
        if isinstance(trjdata,dict):
            self.trjdata = trjdata   # overwritten when load()ing unless none present in saved
        if isinstance(site_properties,numpy.recarray):
            self.site_properties = site_properties  # managed attribute; also builds equivalent_sites_index
        self._filename = filename
        if not (graph is None or properties is None):
            self.graph = NX.DiGraph(name='Transitions between sites')
            self.properties = properties
            self.graph.add_nodes_from(graph.nodes)  # nodes for completeness
            #-------------------------------------------------------------------------------
            # Compute rates for each edge from list of tau (can take a while!);
            # here the format of an edge is defined: (from_site,to_site,rate_tuple)
            ebunch = [(e[0],e[1],self._rate(properties[e]['tau'])) for e in graph.edges]
            self.graph.add_edges_from(ebunch)  # leaves out any nodes that have no edges
            #-------------------------------------------------------------------------------
        elif filename is not None:
            self.load(filename)
        else:
            raise ValueError('provide edges and properties or a filename')

    # convenience functions to access data in an edge (=hop)
    # NX 1.x: edge = (from_site, to_site, {'k': k, 'N':N, 'fit': fit_func})
    # edge from graph.edges(data=True)
    def rate(self,edge):
        """Returns the fastest rate on an edge, in ns^-1"""
        f = self.waitingtime_fit(edge)
        if len(f.parameters) == 1:
            # exp(-k*t)
            k = f.parameters[0]
        elif len(f.parameters) == 3:
            # a exp(-k1*t) + (1-a)exp(-k2*t)
            a,k1,k2 = f.parameters
            ## select rate that contributes more:
            ##if a*k1 > (1-a)*k2:
            if a > 0.5:
                k = k1
            else:
                k = k2
        else:
            raise ValueError('Unknown waitingtime fit, '+str(f))
        return 1000 * k   # 1/ps --> 1/ns

    def number_of_hops(self,edge):
        """Number of transitions recorded."""
        return edge[2]['N']

    def waitingtime_fit(self,edge):
        """Returns the fit function for the edge's waiting time distribution."""
        return edge[2]['fit']

    def from_site(self,edge):
        """Returns the originating site of hop."""
        return edge[0]

    def to_site(self,edge):
        """Returns the site to which a hop is directed."""
        return edge[1]

    def is_from_bulk(self,edge):
        """True if the edge originated in the bulk."""
        return self.from_site(edge) == SITELABEL['bulk']

    def site_properties():
        doc = """Site_properties, indexed by node label.
              Setting this attribut also updates self.equivalent_sites_index."""
        def fget(self):
            return self.__site_properties
        def fset(self,x):
            self.__site_properties = x
            try:
                # recreate equivalent_sites_index (no need to store it)
                Q = x[x.equivalence_label > SITELABEL['bulk']]
                self.equivalent_sites_index = dict(zip(Q.equivalence_label,Q.label))
            except AttributeError:
                # should still work for densities without equivalence sites
                self.equivalent_sites_index = {}
        return locals()
    site_properties = property(**site_properties())

    def compute_site_times(self,verbosity=3):
        """Compute the 'residency' time of each water molecule on each site.

        compute_site_times()

        The 'life time' of a site i is computed as

          theta[i] = <t[*,i]>

        where t[*,i] stands for all waiting times t[j,i] for hops from
        site i to all other sites AND the waiting times t[i] of
        molecules that are not observed to leave site i.

        * The function updates self.theta[site] for each site with an
        array of residency times (in ps).

        * The life times are stored in self.lifetime_avg[site] and
          self.lifetime_std[site]

        TODO:
        * Maybe use the barrier time as well (or a portion thereof,
          perhaps proportional to the barrier height (related to the
          kji) --- rate theory??)
        """
        self.theta = dict()        # dict of arrays
        for edge in self.properties:
            if type(edge) is tuple:
                i,j = edge       # i --> j
            else:
                i = edge         # molecule did not leave site
            # could look at 'tbarrier' as well
            try:
                self.theta[i] = numpy.append(self.theta[i], self.properties[edge]['tau'])
            except KeyError:
                self.theta[i] = self.properties[edge]['tau'].copy()       # make a new array
        # array for accumulated site times
        # assumptions: labels are integers from 0 to smax
        smax = numpy.max(self.graph.nodes())
        self.lifetime_avg = numpy.zeros(smax+1)        # include 0 so that label matches index
        self.lifetime_std = numpy.zeros(smax+1)  # error estimate (biased, N, not N-1)
        for label,thetas in self.theta.items():
            self.lifetime_avg[label] = numpy.average(thetas)
            self.lifetime_std[label] = numpy.std(thetas) # biased

    def compute_site_occupancy(self):
        """Computes occupancies from the residency times theta and updates self.occupancy.

        compute_site_occupancy()

        occupancy::
                        N_i
             o[i] = 1/T Sum theta[i,k]
                        k=1

           where T is the total trajectory time and the sum runs over
           all residency times that were recorded for the site i.

        attributes:

           self.occupancy_avg
                numpy array with occupancies, site label == index
           self.occupancy_std
                numpy array with error estimates for occupancies
                (Delta = Delta(theta)/T; this is a biased estimate because
                Delta(theta) is calculated with N instead of N-1)
        """
        try:
            if len(self.theta) == 0:
                raise ValueError
        except (AttributeError,ValueError):
            # need to compute thetas first
            self.compute_site_times(verbosity=0)

        # array for accumulated site times
        # assumptions: labels are integers from 0 to smax
        smax = numpy.max(self.graph.nodes())
        self.occupancy_avg = numpy.zeros(smax+1)        # include 0 so that label matches index
        self.occupancy_std = numpy.zeros(smax+1)  # error estimate (biased, N, not N-1)
        for label,thetas in self.theta.items():
            self.occupancy_avg[label] = numpy.sum(thetas) / self.trjdata['totaltime']
            self.occupancy_std[label] = numpy.std(thetas) / self.trjdata['totaltime'] # biased

    def save(self,filename=None):
        """Save HoppingGraph as a pickled python object."""
        filename = self.filename(filename,'pickle',set_default=True)
        fh = open(filename,'wb')
        try:
            cPickle.dump(self,fh,cPickle.HIGHEST_PROTOCOL)
        finally:
            fh.close()

    def load(self,filename=None):
        """Reinstantiate HoppingGraph from a pickled HoppingGraph (from save())."""
        filename = self.filename(filename,'pickle',set_default=True)
        fh = open(filename,'rb')
        try:
            h = cPickle.load(fh)
        finally:
            fh.close()
        # restore attributes from the temporary instance (works but not elegant...)
        for attr in 'graph', 'properties', 'filtered_graph', 'trjdata', \
                'node_properties','occupancy_avg','occupancy_std', :
            if hasattr(h,attr):
                self.__dict__[attr] = h.__dict__[attr]
        # managed attributes
        self.site_properties = h.site_properties
        del h

    # XXX: monkey patching, should do this with inheritance/mixin class
    filename = utilities.filename_function

    def filter(self,exclude=None):
        """Create a filtered version of the graph.

        For looking at most things:
        >>> h.filter(exclude={'outliers':True})

        For looking at exchange rates and plotting:
        >>> h.filter(exclude={'outliers':True, 'Nmin':5, 'unconnected':True})

        For export3D do not use the bulk site:
        >>> h.filter(exclude={'outliers':True,'bulk':True})

        This method makes a copy of the hopping graph and applies the
        filter rules to the copy. Other output functions use this copy
        if it exists.

        exclude      dict of components to exclude. May contain
                     {'outliers':True, 'Nmin':integer,
                      'bulk': True, 'unconnected':True}

                     If outliers == True then all edges from the 'outlier' node
                     are deleted previous to displaying the
                     graph. Those edges correspond to particles
                     starting in a region not covered by the intial
                     histogram boundaries and enter a mapped site at a
                     later point in time.

                     With Nmin, any node that has fewer than Nmin transition
                     is discarded.

                     unconnected == True finaly filters all nodes that have no
                     edges left

        """
        G = self.filtered_graph = self.graph.copy()
        if exclude is None:
            return

        try:
            # order matters for filters
            if 'Nmin' in exclude:
                for u,v,p in self.graph.edges(data=True):
                    if p['N'] < exclude['Nmin']:
                        G.remove_edge(u,v)
                G.name += ', N>='+str(exclude['Nmin'])
            if 'outlier' in exclude:
                if exclude['outlier'] and SITELABEL['outlier'] in G:
                    G.remove_node(SITELABEL['outlier'])
            if 'bulk' in exclude:
                if exclude['bulk'] and SITELABEL['bulk'] in G:
                    G.remove_node(SITELABEL['bulk'])
                    G.name += ', bulk site omitted'
            if 'unconnected' in exclude:
                if exclude['unconnected']:
                    delnodes = [n for n in G.nodes if G.degree(n) == 0]
                    G.remove_nodes_from(delnodes)
                    G.name += ', %d sites omitted' % len(delnodes)
                    del delnodes
        except TypeError:
            raise TypeError('"exclude" should have been a dict.')


    def plot_fits(self,ncol=2,nrow=3,dt=None,plottype='log',use_filtered_graph=True,
                  directory='survival_times',format='png',interactive=False,
                  verbosity=3):
        """Plot survival time fit against data.

        plot_fits(ncol=2)

        The time values are taken to cover all measured tau.

        ncol         number of columns
        nrow         number of rows per page
        plottype     'linear' or 'log'
        dt           time step in ps; use value in self.trjdata['dt'] or 1ps
        use_filtered_graph
                     True: use the filtered graph (see filter()),
                     False: use raw data.
        directory    save all pdf files under this directory
        format       file format for plot (png,eps,pdf... depends on matplotlib)
        interactive  False: do not display graphs on scren (default)
                     True: show graphs on screen, can be slow and probably
                     requires ipython as your python shell
        verbosity    chattiness level

        All N graphs are laid out in nrow x ncol grids on as many
        pages/figures as necessary.

        The pages are written as eps/pdf files using a fixed filename
        in the given directory ('survival_times' by default).
        """
        import math
        import os,errno

        import matplotlib
        from .utilities import matplotlib_interactive
        matplotlib_interactive(interactive)
        import matplotlib.pyplot as plt

        set_verbosity(verbosity)

        try:
            os.makedirs(directory)
        except os.error,e:
            if e.errno != errno.EEXIST:
                raise

        # Is there a filtered graph we should use?
        G = self.select_graph(use_filtered_graph)

        if not dt:
            try:
                dt = self.trjdata['dt']
            except (KeyError,AttributeError):
                dt = 1.0    # default guess of dt = 1 ps
        def trange(tau,dt=dt):
            """Linear 1d mesh covering all tau values + 2*dt more"""
            tmin,tmax = 0, numpy.max(tau) + 2*dt
            return numpy.arange(tmin,tmax,dt)
        if plottype is "linear":
            plotfunc = 'plot'
            loc = 'best'
            fntempl = 'S_%02d'
            def tSrange(t_all,S):
                """use all values of t for linear plots"""
                return t_all
        elif plottype is "log":
            plotfunc = 'semilogy'
            loc = 'lower left'
            fntempl = 'S_%02d_log'
            def tSrange(t,S):
                """only use values where S>0"""
                return t[S(t) > 0]
        else:
            raise ValueError('plottype "%s" is not valid' % plottype)

        # extension determines plot output format
        filename_templ = os.path.join(directory,fntempl) + '.' + format

        # Layout all N plots on npage pages in a nrow x ncol grid.
        nplots = G.number_of_edges()
        npage = nrow * ncol
        nfig = int(math.ceil(1.0*nplots/npage))

        if nplots == 0:
            msg = "Filtered graph contains no transitions: no plots available"
            warnings.warn(msg, category=MissingDataWarning)
            logger.warn(msg)
            return

        # erase all figures first
        # (important in interactive sessions with figure pollution)
        for ifig in xrange(1,nfig+1):
            plt.close(ifig)

        font = {'legend': matplotlib.font_manager.FontProperties(size='xx-small'), }
        pm = CustomProgressMeter(nplots, interval=1, offset=1,
                                 format="Preparing page %(other)2d, S(t) plot %(step)3d/%(numsteps)d  [%(percentage)5.1f%%].\r")
        for iplot,e in enumerate(G.edges(data=True)):
            # NX 1.x: e = (n1, n2, data) with data={'k':k, 'N':N, 'fit':fit}
            ifig = iplot/npage + 1
            isub = iplot % npage + 1
            pm.echo(iplot, ifig)

            fig = plt.figure(ifig)  # , figsize=(8,11) letter 8.5" x 11"; use default
            ax = fig.add_subplot(nrow,ncol,isub, frameon=False)

            tau = self.properties[e[0:2]]['tau']
            S = survivalfunction(tau)
            t = trange(tau)
            tS = tSrange(t,S)
            k,N,Sfit = e[2]['k'], e[2]['N'], e[2]['fit']
            if len(Sfit.parameters) == 3:         # get two rates if double exp
                a,k1,k2 = Sfit.parameters
                kns = 1000 * numpy.array([k1,k2]) # rate in 1/ns
                k_legend = "$k_{1}=%.1f/$ns [%.0f%%]\n$k_{2}=%.1f/$ns [%.0f%%]" % \
                           (kns[0],100*a, kns[1],100*(1-a))
            else:
                kns = 1000*k   # rate in 1/ns
                k_legend = r"$k=%.1f/$ns " % kns
            getattr(ax, plotfunc)(tS,S(tS),'k-', tau,S(tau),'k.')
            getattr(ax, plotfunc)(t,Sfit.fit(t), 'r--',linewidth=2)
            ax.legend((r"$S(t)\/ %d\rightarrow%d$" % (e[0],e[1]), r"$N=%d$" % N, k_legend),
                      loc=loc, numpoints=2, prop=font['legend'])
            if ax.is_last_row()  : ax.set_xlabel(r'time $t$/ps')
            if ax.is_first_col() : ax.set_ylabel(r'$S(t)$')

            if ax.is_last_row() and ax.is_last_col():
                filename = filename_templ % ifig
                fig.savefig(filename)
                plt.close(fig)
                # export as eps & convert to pdf externally for older matplotlib
                # (just use the default png!)
                #os.system('epstopdf --outfile=%(pdf)s %(eps)s' % locals())
                #os.remove(eps)
                logger.info("Exported S(t) plot page %(ifig)d to file %(filename)s.", locals())

    def tabulate_k(self):
        """List of tuples (from, to, rate (in 1/ns), number of transitions)."""
        return [(self.from_site(edge), self.to_site(edge),
                 self.rate(edge), self.number_of_hops(edge) ) for edge in self.graph.edges(data=True)]

    def show_rates(self, filename=None):
        """Print the rates (in 1/ns) between sites, and the total number of observations.

        show_rates(file=filename)

        By default, prints to stdout but if *file* = filename then
        filename is opened and data are written to the file.

        A description of the fit function used to obtain the rate is
        also printed in the last column.

        Only the "dominant" rate is shown; see the fit_func
        description for cases when two rates were computed.

        .. SeeAlso:: :func:`HoppingGraph.tabulate_k`.
        """
        if not filename is None:
            stream = open(filename, 'w')
        else:
            stream = sys.stdout

        try:
            for edge in self.graph.edges(data=True):
                line = "%3d --> %3d   %6.1f  1/ns  %3d  %s\n" % \
                    (self.from_site(edge), self.to_site(edge), self.rate(edge),
                     self.number_of_hops(edge), str(self.waitingtime_fit(edge)))
                stream.write(line)
        finally:
            if not filename is None:
                stream.close()

    def equivalent_sites_stats(self,elabels,equivalence=True):
        """Statistics about one or a list of equivalence sites.

        g.equivalent_sites_stats(elabels,equivalence=True)

        :Arguments:
        elabels         a single label or a list of node labels
        equivalence     True: interprete elabels as 'equivalence labels', i.e. the label
                        attached to a site common to two densities
                        False: elabels are labels local to the graph
        """
        # TODO: return as a dict for processing
        #       add graph attribute (default to 0 or None - or a reference to self?)
        try:
            elabels[0]
        except TypeError:
            elabels = [elabels]
        if equivalence:
            nbunch = [self.equivalent_sites_index[l] for l in elabels]
        else:
            nbunch = elabels
        for node in nbunch:
            print "%s   elabel=%d   node=%d" % \
                (self.site_properties.equivalence_name[node],
                 self.site_properties.equivalence_label[node],
                 node)
            print "Node: "
            try:
                print self.properties[node]
            except KeyError:
                print "---"
            print "In-edges:"
            print self.graph.in_edges(node)
            print "Out-edges:"
            print self.graph.out_edges(node)
            print "\n"

    def stats(self,data=None):
        """Statistics for the hopping graph.

          stats([data=dict]) --> dict

        Without the data argument, the method just returns some
        interesting values gathered from the graph and the density. If
        a data dictionary is given, then the raw data are loaded into
        the dict and can be processed further by histogramming etc.

        :Arguments:
           data
               optional dictionary to hold raw data for processing;
               modified by method

        :Returns:  dictionary with expressive keys, holding the results
        """
        if (not hasattr(self,'site_properties') or self.site_properties is None):
            raise AttributeError('Stats require site_properties annotation.')
        graph = self.graph   # makes little sense to use filtered graph
        stats = {}

        # get stats for all nodes; note that compound equivalence
        # sites in graph have replaced original sites
        nodes = graph.nodes()
        nodes.remove(SITELABEL['bulk'])  # exclude bulk from stats
        nodes = numpy.asarray(nodes)

        #------------------------------------------------------------
        # Graph stats
        #------------------------------------------------------------
        # general stats including bulk
        stats['G_order'] =  graph.order()
        stats['G_edges'] = graph.number_of_edges()
        stats['G_degree'] = 2.0*graph.number_of_edges()/graph.order()
        stats['G_degree_in'] = stats['G_degree_out'] = stats['G_degree']/2.0
        stats['G_degree_max'] = numpy.max(graph.degree().values())
        stats['G_degree_min'] = numpy.min(graph.degree().values())
        stats['G_degree_in_max'] = numpy.max(graph.in_degree().values())
        stats['G_degree_in_min'] = numpy.min(graph.in_degree().values())
        stats['G_degree_out_max'] = numpy.max(graph.out_degree().values())
        stats['G_degree_out_min'] = numpy.min(graph.out_degree().values())

        # general stats excluding bulk and bulk --> site edges
        stats['G_order_nobulk'] =  graph.order() - 1
        stats['G_edges_nobulk'] = stats['G_edges'] - 0.5*graph.degree(SITELABEL['bulk'])
        stats['G_degree_nobulk'] = 2.0*stats['G_edges_nobulk']/stats['G_order_nobulk']
        stats['G_degree_nobulk_max'] = numpy.max(graph.degree(nodes).values())
        stats['G_degree_nobulk_min'] = numpy.min(graph.degree(nodes).values())
        stats['G_degree_in_nobulk'] = numpy.average(graph.in_degree(nodes).values())
        stats['G_degree_out_nobulk'] = numpy.average(graph.out_degree(nodes).values())
        stats['G_degree_in_nobulk_max'] = numpy.max(graph.in_degree(nodes).values())
        stats['G_degree_out_nobulk_max'] = numpy.max(graph.out_degree(nodes).values())
        stats['G_degree_in_nobulk_min'] = numpy.min(graph.in_degree(nodes).values())
        stats['G_degree_out_nobulk_min'] = numpy.min(graph.out_degree(nodes).values())

        internal = self.internal_sites()
        stats['G_internal'] = len(internal)     # includes isolated sites
        isolated = self.isolated_sites()
        stats['G_isolated'] = len(isolated)

        try:
            data['G_nodes'] = nodes
            data['G_degree'] = graph.degree(nodes).values()         # NX 1.x: returns a dict, to
            data['G_degree_in'] = graph.in_degree(nodes).values()   # keep consistent with old behaviour of hop
            data['G_degree_out'] = graph.out_degree(nodes).values() # add flat lists to data with values()
            data['G_internal'] = internal
            data['G_isolated'] = isolated
        except TypeError:
            pass

        #------------------------------------------------------------
        # kinetics stats
        #------------------------------------------------------------
        if not hasattr(self,'lifetime_avg'):
            self.compute_site_times()
        if not hasattr(self,'occupancy_avg'):
            self.compute_site_occupancy()
        # lifetime (!= residence time)
        stats['site_lifetime_avg'] = numpy.average(self.lifetime_avg[nodes])
        stats['site_lifetime_med'] = numpy.median(self.lifetime_avg[nodes])
        stats['site_lifetime_std'] = numpy.std(self.lifetime_avg[nodes])
        stats['site_lifetime_max']= numpy.max(self.lifetime_avg[nodes])
        stats['site_lifetime_min']= numpy.min(self.lifetime_avg[nodes])
        # occupancy from life times
        stats['site_occupancy_kin_avg'] = numpy.average(self.occupancy_avg[nodes])
        stats['site_occupancy_kin_med'] = numpy.median(self.occupancy_avg[nodes])
        stats['site_occupancy_kin_std'] = numpy.std(self.occupancy_avg[nodes])
        stats['site_occupancy_kin_max'] = numpy.max(self.occupancy_avg[nodes])
        stats['site_occupancy_kin_min'] = numpy.min(self.occupancy_avg[nodes])

        try:
            data['site_lifetime_avg'] = self.lifetime_avg[nodes]
            data['site_lifetime_std'] = self.lifetime_std[nodes]
            data['site_occupancy_kin_avg'] = self.occupancy_avg[nodes]
            data['site_occupancy_kin_std'] = self.occupancy_std[nodes]
        except TypeError:
            pass

        #------------------------------------------------------------
        # density stats
        #------------------------------------------------------------
        sp = self.site_properties[nodes]
        stats['site_volume_avg'] = numpy.average(sp.volume)
        stats['site_volume_med'] = numpy.median(sp.volume)
        stats['site_volume_avg'] = numpy.std(sp.volume)
        stats['site_occupancy_rho_avg'] = numpy.average(sp.occupancy_avg)
        stats['site_occupancy_rho_med'] = numpy.median(sp.occupancy_avg)
        stats['site_occupancy_rho_std'] = numpy.std(sp.occupancy_avg)

        # re-derive N_equivalence_sites and N_subsites from site_properties
        # see hop.sitemap.stats
        try:
            stats['site_N_equivalence_sites'] = len(sp[sp.equivalence_label != 0])
            stats['site_N_subsites'] = len(self.site_properties[self.site_properties.equivalence_site != 0])
        except AttributeError:
            # site_properties not fully defined; should perhaps raise...
            warnings.warn('site_properties not complete: %s misses equivalence site data' %
                          self.filename(),       # let's hope filename is set...
                          category=MissingDataWarning)
            stats['site_N_equivalence_sites'] = 0
            stats['site_N_subsites'] = 0

        try:
            data['site_volume'] = sp.volume
            data['site_occupancy_rho_avg'] = sp.occupancy_avg
        except TypeError:
            pass

        #------------------------------------------------------------
        # rates stats
        #------------------------------------------------------------
        # * collect all rates except bulk --> x
        # * sort fast and slow rates separately (2 exp) and single rates (exp)
        # TODO

        return stats

    def internal_sites(self):
        """Returns list of sites that have no connection to the bulk."""
        nodes = self.graph.nodes()
        nodes.remove(SITELABEL['bulk'])  # exclude bulk from stats, not strictly necessary
        nodes = numpy.asarray(nodes)
        is_internal = numpy.array([self.is_internal(n) for n in nodes])
        return asiterable(nodes[is_internal])

    def isolated_sites(self):
        """Returns list of sites that have no other connections."""
        nodes = numpy.asarray(self.graph.nodes())
        is_isolated = (numpy.array(self.graph.degree(nodes)) == 0)
        return asiterable(nodes[is_isolated])

    def is_connected(self,n1,n2):
        """True if node n1 has any connection to the site n2."""
        return self.graph.has_successor(n1,n2) or self.graph.has_predecessor(n2,n1)

    def is_internal(self,n):
        """True if site n has no connection to the bulk."""
        return not self.is_connected(n,SITELABEL['bulk'])

    def is_isolated(self,n):
        """True if site n has no connections to other sites (ie its degree equals 0)."""
        return self.graph.degree(n) == 0

    def connectedness(self,n):
        """Return values that measure connectedness (can be used in occupancy field)"""
        # XXX: quick hack... should be table-driven
        if self.is_isolated(n):
            return 0.01
        if self.is_internal(n):
            return 0.5
        return 1.0

    def rates(self,n,use_filtered_graph=True):
        """Returns k_tot, k_in, k_out (and N_*) for site n (bulk rates omitted from k).

        dictionary = rates(n,use_filtered=True)

        k_in  = sum_j k_nj  (> 0)     j<>bulk
        k_out = sum_j k_jn  (< 0)     j<>bulk
        k_tot = k_in + k_out

        Note that k_tot should be ~0 if a bulk rate is included
        because the graph should obey detailed balance.
        """
        G = self.select_graph(use_filtered_graph=use_filtered_graph)
        k_in  = numpy.sum( [ self.rate(e) for e in G.in_edges(n, data=True)
                             if not self.is_from_bulk(e)] )
        k_out = numpy.sum( [-self.rate(e) for e in G.out_edges(n, data=True)] )
        k_tot = k_in + k_out
        N_in  = numpy.sum( [ self.number_of_hops(e) for e in G.in_edges(n, data=True)] )
        N_out = numpy.sum( [-self.number_of_hops(e) for e in G.out_edges(n, data=True)] )
        N_tot = N_in + N_out
        ratedict={'k_tot':k_tot,'k_in':k_in,'k_out':k_out,
                'N_tot':N_tot,'N_in':N_in,'N_out':N_out}
        return  ratedict

    def show_site(self,sites,use_filtered_graph=True):
        """Display data about sites (list of site labels or single site)."""
        if not iterable(sites):
            sites = [sites]
        for s in sites:
            self._show_site(s,use_filtered_graph=use_filtered_graph)
            msg(0,'\n')

    def _show_site(self,site,use_filtered_graph=True):
        """Display data about a single site."""
        G = self.select_graph(use_filtered_graph)
        rule = "-"*50+"\n"
        msg(0,"Statistics for site %d" % site)
        equivalence_name = self.site_properties.equivalence_name[site].strip()
        if equivalence_name:
            msg(0," (common site '%s')" % equivalence_name)
        msg(0,"\n")
        if use_filtered_graph:
            msg(0,"(Filtered graph '%s')\n" % G.name)
        sp = self.site_properties[site]
        msg(0,"  occupancy    = %g\n" % sp.occupancy_avg)
        msg(0,"  site volume  = %g A^3\n" % sp.volume)
        msg(0,"  internal site = %s    isolated site = %s\n" % (self.is_internal(site),
                                                                self.is_isolated(site)))

        msg(0,rule)
        msg(0,"%11s %12s  1/ns  %3s\n" % ("hop","rate in","N"))
        msg(0,rule)
        msg(0,"In-edges:\n")
        for e in G.in_edges(site, data=True):
            self.show_edge(e)
        msg(0,"Out-edges:\n")
        for e in G.out_edges(site, data=True):
            self.show_edge(e)
        msg(0,"Total rates, excluding bulk rate k(%d-->%d):\n" % \
                (SITELABEL['bulk'],site))
        rates = self.rates(site,use_filtered_graph)
        fmt = "%11s %+12.1f  1/ns  %+4d\n"
        msg(0,fmt % ('k_in',rates['k_in'],rates['N_in']))
        msg(0,fmt % ('k_out',rates['k_out'],rates['N_out']))
        msg(0,rule)
        msg(0,fmt % ('k_tot',rates['k_tot'],rates['N_tot']))
        msg(0,rule)

    def show_edge(self,edge):
        msg(0,"%3d --> %3d %12.1f  1/ns  %3d\n" %
            (self.from_site(edge), self.to_site(edge), self.rate(edge),
             self.number_of_hops(edge) ))

    def show_total_rates(self,use_filtered_graph=True):
        """Display total rates for all nodes (excluding bulk --> site contributions)."""
        msg(0,"Total flux in ns^-1 for each site (excluding bulk --> site)\n")
        if use_filtered_graph:
            msg(0,"(Filtered graph '%s')\n" % self.select_graph(use_filtered_graph).name)
        for site in self.graph:
            if site == SITELABEL['bulk']: continue
            rates = self.rates(site,use_filtered_graph)
            if rates['k_tot'] > 0:   marker = '*'
            else:           marker = ' '
            rates.update({'site':site,'marker':marker})
            msg(0,"%(site)2d  %(marker)s k_tot=%(k_tot)+9.1f    k_in=%(k_in)9.1f   k_out=%(k_out)9.1f\n" % rates)


    def export(self,filename=None,format='XGMML',use_filtered_graph=True,
                    use_mapped_labels=True):
        """Export the graph to a graph format or an image.

        export('hopgraph',format='XGMML',use_filtered_graph=True)

        :Arguments:

        filename     name for the output files; appropriate suffixes are added
                     automatically
        format       output format: graphs (XGMML or DOT) or image  (png, jpg, ps, svg)
        use_filtered_graph
                     By default, the filtered graph (see the filter() method) is
                     plotted. If set to False then the original HoppingGraph is
                     used instead.
        use_mapped_labels
                     If site_properties is provided then each node that has been
                     identified to exist in a reference network is coloured black
                     and the mapped label is printed instead of the graph label.
        """
        __exporters = {'XGMML': self._export_xgmml,
                       'DOT': self._export_dot,
                       'ps': self._export_image,
                       'svg': self._export_image,
                       'png': self._export_image,
                       'jpg': self._export_image,
                       }
        try:
            __exporters[format](filename=filename,use_filtered_graph=use_filtered_graph,
                                format=format,use_mapped_labels=use_mapped_labels)
        except KeyError:
            raise ValueError("Format "+str(format)+" not supported, choose one of "+
                             str(__exporters.keys()))


    def _export_dot(self,filename=None,use_filtered_graph=True,use_mapped_labels=True,
                    **otherargs):
        """Export the graph as dot file.

        export('graph')

        :Arguments:

        filename     name for the output files; appropriate suffixes are added
                     automatically
        use_filtered_graph
                     By default, the filtered graph (see the filter() method) is
                     plotted. If set to False then the original HoppingGraph is
                     used instead.
        use_mapped_labels
                     If site_properties is provided then each node that has been
                     identified to exist in a reference network is coloured black
                     and the mapped label is printed instead of the graph label.

        The graph is only written to an image file if an image format
        is supplied. See
        https://networkx.lanl.gov/reference/pygraphviz/pygraphviz.agraph.AGraph-class.html#draw
        for possible output formats but png, jpg, ps are safe bets.

        See http://graphviz.org/doc/info/attrs.html for attributes.

        Note: On Mac OS X 10.3.9+fink the pygraphviz rendering is
        buggy and does not include node labels. Simply use the
        exported .dot file and use Mac OS X graphviz from
        http://www.pixelglow.com/graphviz/
        """
        dotname = self.filename(filename,'dot')       # load into graphviz
        G = self._make_AGraph(use_filtered_graph=use_filtered_graph,use_mapped_labels=use_mapped_labels)
        G.write(dotname)

    def _export_image(self,filename=None,format='png',use_filtered_graph=True,use_mapped_labels=True,
                      **otherargs):
        G = self._make_AGraph(use_filtered_graph=use_filtered_graph,use_mapped_labels=use_mapped_labels)
        gfxname = self.filename(filename,format)
        G.layout(prog=otherargs.pop('layouter', 'dot'))
        G.draw(path=gfxname,format=format)  # buggy on Mac OS X

    def _make_AGraph(self,use_filtered_graph=True,use_mapped_labels=True):
        """Returns AGraph ('attributed graph') of the filtered or current graph."""

        import pygraphviz

        if use_mapped_labels and \
           (not hasattr(self,'site_properties') or self.site_properties is None):
            raise AttributeError('Mapped site labels require site_properties annotation.')

        graph = self.select_graph(use_filtered_graph)

        G = pygraphviz.AGraph(directed=True)
        G.add_nodes_from(graph.nodes())
        G.add_edges_from(graph.edges) # n1 n2

        # attributes: http://www.graphviz.org/doc/info/attrs.html
        # colors: http://www.graphviz.org/doc/info/colors.html
        G.graph_attr['label'] = graph.name + ', k in 1/ns (N)'
        G.graph_attr['splines'] = 'true'
        G.graph_attr['fontcolor'] = 'black'
        G.node_attr['fontcolor'] = 'white'
        G.node_attr['color'] = 'black'
        G.node_attr['shape'] = 'circle'
        G.node_attr['fillcolor'] = 'gray76'    # middle gray
        G.node_attr['style'] = 'filled'
        G.edge_attr['color'] = 'gray36'        # darker gray
        G.edge_attr['dir'] = 'forward'

        # attach k_ji to edge (printed in 1/ns), number of transitions in ()
        for e in G.edges(data=True):    # TODO: NX 1.x : do we need data=True here?
            u,v = map(int,e[:2])        # need nodes as integers
            k,N,fit = e[2]['k'], e[2]['N'], e[2]['fit']   # extract k from HoppingGraph
            if len(fit.parameters) == 3:         # get two rates if double exp
                a,k1,k2 = fit.parameters
                kns = 1000 * numpy.sort([k1,k2]) # rate in 1/ns
                k_dominant = kns[1]              # only display fastest rate
            else:
                k_dominant = 1000*k              # rate in 1/ns
            k_legend = "%(k_dominant).1f (%(N)d)" % locals()
            e.attr['headlabel'] = k_legend
            if use_mapped_labels:
                for k in u,v:
                    node = G.get_node(k)
                    label = int(node)
                    try:
                        mapped_label = self.site_properties.equivalence_name[label]
                        if mapped_label.strip():
                            node.attr['label'] = mapped_label.strip()
                            node.attr['fillcolor'] = 'red'
                            node.attr['fontcolor'] = 'white'
                    except ValueError:
                        pass
        return G


    def _export_xgmml(self,filename=None,use_filtered_graph=True,**otherargs):
        # quick hack
        G = self.select_graph(use_filtered_graph)
        filename = self.filename(filename,'xgmml')
        xattr = {'id':'HoppingGraph', 'label':G.name, 'description':G.name,
                 'softwareVersion': hop.__version__,
                 }
        xml = open(filename,'w')
        xml.write("""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n"""
                  """<!DOCTYPE graph PUBLIC "-//John Punin//DTD graph description//EN" "http://www.cs.rpi.edu/research/groups/pb/punin/public_html/XGMML/GML_XGMML/xgmml.dtd">\n"""
                  """<graph label="%(label)s" directed="1" id="%(id)s" """
                  """xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" """
                  """xmlns:dc="http://purl.org/dc/elements/1.1/" """
                  """xmlns:xlink="http://www.w3.org/1999/xlink" """
                  """xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" """
                  """xmlns="http://www.cs.rpi.edu/XGMML" """
                  """Vendor="hop.graph.HoppingGraph._export_xgmml">\n"""
                  """\t<att name="documentVersion" value="1.0"/>\n"""
                  """\t<att name="softwareVersion" value="%(softwareVersion)s"/>\n"""
                  """\t<att name="description" value="=%(description)s"/>\n""" % xattr)

        # distance of site center from geometrical center
        # (a bit hackish so that we only include real sites)
        sitecenters = numpy.array([v for v in
                                   self.site_properties.center[SITELABEL['bulk']+1:]])
        center = numpy.mean(sitecenters,axis=0)
        cdist = numpy.sqrt(numpy.sum((sitecenters - center)**2,axis=1))
        # add back interstitial and bulk so that site label directly indexes into cdist
        centerdistance = numpy.concatenate( ([0.0]*(SITELABEL['bulk']+1),cdist) )
        # del cdist, center, sitecenters
        for k in G.nodes():
            label = str(k)
            site = k
            xattr = {'id':site,'label':label,
                     'equivalence_label':self.site_properties.equivalence_name[site],
                     'volume':self.site_properties.volume[site],
                     'occupancy_avg':self.site_properties.occupancy_avg[site],
                     'distance':centerdistance[site],
                     'has_bulkconnection':self.is_connected(site,SITELABEL['bulk']),
                     'rates':self.rates(site)['N_tot'],
                        }
            if xattr['equivalence_label']:
                xml.write("""\t<node id="%(id)d" label="%(label)s/%(equivalence_label)s">\n""" % xattr)
            else:
                xml.write("""\t<node id="%(id)d" label="%(label)s">\n""" % xattr)
            xml.write("""\t\t<att type="real" name="volume" value="%(volume)g"/>\n""" % xattr)
            xml.write("""\t\t<att type="real" name="occupancy" value="%(occupancy_avg)g"/>\n""" % xattr)
            xml.write("""\t\t<att type="real" name="distance" value="%(distance)g"/>\n"""  % xattr)
            xml.write("""\t\t<att type="integer" name="has_bulkconnection" value="%(has_bulkconnection)d"/>\n""" % xattr)

            xml.write("""\t\t<att type="real" name="rates" value="%(rates)r"/>\n""" % xattr)
            ### xml.write("""\t\t<att type="" name="" value=""/>\n""")
            xml.write("""\t</node>\n""")
        for e in G.edges(data=True):
            u,v = self.from_site(e), self.to_site(e)   # === u,v = e[:2]
            xattr = {'u':u,'v':v,'label':"%d -> %d" % (u,v),'id':"%d -> %d" % (u,v)}
            xattr['rate'] = self.rate(e)
            xattr['N'] = self.number_of_hops(e)
            xml.write("""\t<edge source="%(u)d" target="%(v)d" label="%(label)s">\n""" % xattr)
##                      """label="%(label)s" id="%(id)s">\n""" % xattr)
            xml.write("""\t\t<att type="string" name="EDGE_TYPE" value="DirectedEdge"/>\n""")
            xml.write("""\t\t<att type="string" name="interaction" value="DirectedEdge"/>\n""")
            xml.write("""\t\t<att type="integer" name="N" value="%(N)d"/>\n""" % xattr)
            xml.write("""\t\t<att type="real" name="rate" value="%(rate)f"/>\n""" % xattr)
            xml.write("""\t</edge>\n""")
        xml.write("""</graph>\n""")
        xml.close()


    def export3D(self,density=None,filename=None,use_filtered_graph=True):
        """Export pdb and psf file for visualization in 3D.

        >>> h.export3D()
        Uses h.site_properties if it exists.

        >>> h.export3D(density)
        Uses a (hopefully matching) Density object to pull in site_properties.


        :Arguments:
        density      hop.sitemap.Density with full site_properties
        filename     prefix for output files: <filename>.psf and <filename>.pdb
        use_filtered_graph
                     define a filtered graph with h.filter() first

        The method writes a psf and a pdb file from the graph, suitable
        for visualization in, for instance, VMD.

        Sites are represented as residues of resname 'NOD'; each site
        is marked by one 'ATOM' (of type CA) at the center of geometry
        of the site. Edges are bonds between those pseudo atoms.

        #Currently: B-factor 1 if common site label exist, 0 otherwis
        #           occupancy: avg site occupancy
        # (but this should become customizable)

        One should use a filtered graph with the bulk site removed for
        visualization.

        Bugs:
        * with a filtered graph, the degree is the one of the filtered
          graph and not of the real underlying graph
        * cannot yet select what to display in B-factor and occupancy field:
          choose from: ['identity','occupancy','degree','volume']
        """
        if density is None:
            if not hasattr(self,'site_properties') or self.site_properties is None:
                raise AttributeError('Requires site_properties or a Density instance.')
            props = self.site_properties
        else:
            try:
                props = density.site_properties
            except AttributeError:
                raise TypeError('density must be a hop.sitemap.Density instance')
        graph = self.select_graph(use_filtered_graph)

        self.write_psf(graph,props,filename)        # atoms are numbered consecutively...
        self.write_pdb(graph,props,filename)  # ..and residues correspond to sites

    def write_pdb(self,graph,props,filename=None):
        # TODO: replace with MDAnalysis PDB writer

        B = Bio.PDB.StructureBuilder.StructureBuilder()
        B.init_structure('graph')
        B.init_model(0)
        B.init_seg('GRPH')
        B.init_chain('A')
        # note that empty fields MUST be written as a blank ' ' (ie a single space)
        # or they are interpreted as altLoc specifiers named '' (weird...)
        # ATOM numbering is done consecutively here and in write_psf() and
        # there is only a single atom per residue.
        for node in graph:   # node is the label==resid and it must be an integer
            pos = props[node].center
            vol = props[node].volume
            occ = props[node].occupancy_avg
            degree = graph.degree(node)           ## TODO w/filtered (may be off by 1)
            # write common atoms as N, different as CA (HACK...)
            commonlabel = props.equivalence_name[node].strip()
            if commonlabel:
                identity = 1.0
                aname = 'N'
                atype = 'N'
            else:
                identity = 0.0
                aname = 'CA'       # choose the same identifiers as in pdb
                atype = 'CA'       # choose the same identifiers as in pdb
            pdb_occupancy = occ     # should be able to choose from volume, occupancy, degree, identity
            # connectedness: 1: standard site, 0.5: internal (no bulk), 0.01: isolated
            # (should do this as separate residues)
            pdb_beta = self.connectedness(node)
            B.init_residue('NOD',' ',node,' ') # choose same identifiers as in write_psf
            try:
                B.init_atom(aname,pos,pdb_beta,pdb_occupancy,' ',atype, element='C')
            except TypeError:
                # newer version of Biopython ?
                B.init_atom(aname,pos,pdb_beta,pdb_occupancy,' ',atype)
        io=Bio.PDB.PDBIO()
        s = B.get_structure()
        io.set_structure(s)
        pdbfile = self.filename(filename,'pdb')
        io.save(pdbfile)

    def write_psf(self,graph,props,filename=None):
        """Pseudo psf with nodes as atoms and edges as bonds"""
        # Standard no CHEQ format for a Charmm PSF file:
        psf_ATOM_format = '%(iatom)8d %(segid)4s %(resid)-4d %(resname)4s '+\
                          '%(name)-4s %(type)4s %(charge)-14.6f%(mass)-14.4f%(imove)8d\n'

        psffilename = self.filename(filename,'psf')
        psf = open(psffilename,'w')
        psf.write('PSF\n\n')
        psf.write('%7d !NTITLE\n' % 2)
        psf.write('* Graph topology written by\n'+\
                  '* hop {0}\n'.format(hop.__version__))
        psf.write('\n')

        # ATOMS
        psf.write('%6d !NATOM\n' % graph.number_of_nodes())
        segid = 'GRPH'     # choose the same identifiers as in pdb
        resname = 'NOD'    # choose the same identifiers as in pdb
        charge = 0
        mass = 1.0
        imove = 0            # no fixed 'atoms'
        node2iatom = dict()  # needed for BONDS
        for iatom,node in enumerate(graph.nodes):
            # atom numbering starts at 1, so iatom+1
            iatom += 1
            node2iatom[node] = iatom
            # write common atoms as N, different as CA (HACK...)
            commonlabel = props.equivalence_name[node].strip()
            if commonlabel:
                identity = 1.0
                aname = 'N'
                atype = 'N'
            else:
                identity = 0.0
                aname = 'CA'       # choose the same identifiers as in pdb
                atype = 'CA'       # choose the same identifiers as in pdb
            psf.write(psf_ATOM_format %
                      {'iatom':iatom, 'segid':segid, 'resid':node,
                       'resname':resname, 'name':aname, 'type':atype,
                       'charge':charge, 'mass':mass,'imove':imove} )
        # BONDS: fortran fmt03='(8I8)'
        # (note: write directed bonds and one atom per residue/node)
        psf.write('%6d !NBOND: bonds\n' % graph.number_of_edges())
        for n,(i,j,p) in enumerate(graph.edges(data=True)):
            psf.write('%8i%8i' % (node2iatom[i],node2iatom[j]))
            if (n+1) % 4 == 0: psf.write('\n')

        # ignore all the other sections (don't make sense anyway)
        psf.close()

    def select_graph(self,use_filtered_graph):
        """Returns filtered graph for True argument, or the raw graph otherwise)"""
        if use_filtered_graph:
            # Is there a filtered graph we should be using?
            try:
                graph = self.filtered_graph
                graph.number_of_nodes()
            except AttributeError:
                raise ValueError('No filtered graph defined; create one with %s.filter().' %
                                 self.__class__.__name__)
        else:
            graph = self.graph
        return graph

    def _rate(self,taus,method='survivalfunction',block_w=200):
        """Compute the rate i--> j, k_ji, from the hopping times tau_ji.

        :Returns: dict(kij=*kji*, N=*N*, fit=*fit*)

                  - *kji* : rate constant in 1/ps
                  - *N* : number of events
                  - *fit* : instance of :class:`fit_func`

        :Arguments:
            *taus*
                array of all waiting times
            *method*
               'survivalfunction'
                    compute as a fit of a*exp(-k1*t)+(1-a)*exp(-k2*t) or exp(-k*t)
                    to the survival function S(t); the double exp is tried first
                    and then decided heuristically if a single exp is a better choice.
                    Heuristics: Use single exp if

                    * number of data points is <= 10
                    * double exp asymmetry abs(0.5-a) > 0.49
                    * k1<0 or k2<0

        :Bugs:
          - Does not switch to single exponential if double exponential fit fails
            to converge.

        .. Notes:: Should probably use the integral of the double-exponential fit as an
          approximation for the rate constant instead of just using the slower
          one (G. Hummer, pers. comm.)

        :TODO: Implement k from integral as new method 'integral' and also
               compare to the calculation from the exponential fits::

                single exp:
                        f(t) = exp(-kt)
                        k = [integral_0^\infty dt f(t)] ^ -1 = k

                double exponential:
                        f(t) = a exp(-k1 t) + (1-a) exp(-k2 t)
                        k_comb =  [integral_0^\infty dt f(t)] ^ -1 = [a/k1 + (1-a)/k2]^-1

        """

        if method is 'survivalfunction':
            # TODO: use time step for sampling step ?
            # TODO: coarser dt for long runs, eg dt = 10 ... 100!!
            #       This can blow up in memory in survivalfunction() because of too many t (len(x) big!)
            #       or use an interpolated survival function
            dt = 1.0                                 # 1.0 == dt in sim  # TODO: adapt dt
            tmax = numpy.max(taus) + 1.0*dt
            N = len(taus)                            # perhaps sample N times from S ???
            x = numpy.arange(0,tmax+1,dt)            # ... just do lin mesh
            S = survivalfunction(taus,block_w=block_w)
            enough_data = (len(taus) > 10)
            is_single_exp = False                    # try double exp first
            if enough_data:                          # only attempt exp2 with >10 data points
                fit = fitExp2(x,S(x))                # a*exp(-k1*x) + (1-a)*exp(-k2*x)
                a,k1,k2 = fit.parameters
                # heuristics to decide if double exp fit is good
                is_single_exp = abs(0.5 - a) > 0.49  #=0.5 - 0.01  # very asymmetric fit
                all_k = numpy.array([k1,k2])
                if numpy.any(all_k < 0):             # at least one k<0: no sum of exp
                    is_single_exp = True
                else:
                    k = all_k[all_k > 0].max()       # select dominant (fast) rate
            if not enough_data or is_single_exp:
                fit = fitExp(x,S(x))                 # exp(-k*x)
                k, = fit.parameters
        else:
            raise NotImplemented('Method "%s" is not implemented.' % method)
        return dict(k=k, N=N, fit=fit)

class CombinedGraph(HoppingGraph):
    """Hybrid graph between hop graphs that share common nodes."""
    _export_formats = ['XGMML','dot']     # graph formats in export()

    def __init__(self,g0=None,g1=None,filename=None):
        if not (g0 is None or g1 is None):
            if not (isinstance(g0,HoppingGraph) and isinstance(g1,HoppingGraph)):
                raise TypeError("g0 and g1 must be <hop.graph.HoppingGraph>s.")
            # build combined graph
            # 1) add nodes to new graph. Use site_properties.equivalence_name or '1.LABEL', '2.LABEL'
            #    --> make a dict for label --> combined graph label
            # 2) add edges (allow multiple edges) and mark up edge
            #    (can I do multiple parallel edges with pygraphviz ??)
            self._filename = filename
            self.graph = NX.MultiDiGraph(name='Hybrid graph')
            self.graphs = [g0,g1]  # all graphs
            self.g_labels = None   # shape = (2,size of combined graph), one col per graph
                          # g_labels[igraph,combined_label] --> original label in igraph
            self.node_properties = None # will become a numpy.recarray
            self.edge_properties = {}   # indexed by (u,v) = dict(common=common?, plist=[p0, p1])

            # fix: bulk sites should be common, too:
            self.graphs[0].site_properties.equivalence_name[SITELABEL['bulk']] = 'bulk'
            self.graphs[1].site_properties.equivalence_name[SITELABEL['bulk']] = 'bulk'

            # fix: in rare circumstances, an outlier (-1) hop (1,-1) or (-1,1) is recorded
            # so we prune the outlier node --- TODO: check trajectory _coord2hop !!
            for g in self.graphs:
                try:
                    g.graph.remove_node(SITELABEL['outlier'])
                except NX.NetworkXError:
                    pass

            # Set up temporary list to build arrays, start with interstitial at pos = 0
            # tmp list for graph labels with (g0_label, g1_label) for each new label:
            g_labels = [numpy.array([None]*len(self.graphs))]
            # temp list to build recarray: new label, descr, common?, graphmask
            properties = [(SITELABEL['interstitial'],'interstitial',False,0)]
            ICOMLABL = 1               # index for common label in properties at pos 1
            ICOMMON = 2                # True if common to all graphs, False otherwise
            IGRAPHS = 3                # bit field (mask) containing graph numbers as sum 2**igraph
            commonlabel2newnode = dict()  # helper for finding common sites
            inode = SITELABEL['bulk']     # new node labels for hybrid graph start at 'bulk'
            for ig,g in enumerate(self.graphs):
                node2commonlabel = dict(zip(g.site_properties.label,
                                            g.site_properties.equivalence_name))
                for node in g.graph:
                    commonlabel = node2commonlabel[node].strip()
                    newnode = inode
                    inode += 1
                    if commonlabel:
                        if commonlabel not in commonlabel2newnode:
                            self.graph.add_node(newnode)
                            commonlabel2newnode[commonlabel] = newnode  # remember common site's inode
                        else:
                           newnode = commonlabel2newnode[commonlabel]
                           inode -= 1  # after all, don't need a new node
                    try:
                        properties[newnode] # build properties[] in parallel with g_labels[]
                    except IndexError:
                        properties.append([newnode,None,False,0])   # basic record: node,com.label,common?,g
                        g_labels.append(numpy.array([None]*len(self.graphs)))
                    g_labels[newnode][ig] = node  # lookup for node in graph ig
                    properties[newnode][ICOMLABL] = commonlabel
                    properties[newnode][ICOMMON] = commonlabel != '' # convenience flag
                    properties[newnode][IGRAPHS] |= 2**ig            # bit field indicating the original graph
            self.node_properties = numpy.rec.fromrecords(properties,names='label,equivalence_name,common,graphmask')
            self.g_labels = numpy.array(g_labels).transpose()
            del properties
            del g_labels

            site_properties = []  # temporary list to build rec array
            for np in self.node_properties:
                # all but interstitial == [None,None] should have at least one entry in gNnodes
                gNnodes = dict([(ig,g_label) \
                                for ig,g_label in enumerate(self.g_labels[:,np.label]) if g_label])
                volume = 0
                occupancy_avg = 0
                occupancy_std = 0
                center = numpy.array([0.0,0.0,0.0])
                distance = 0.0
                for ig,gNnode in gNnodes.items():
                    sp = self.graphs[ig].site_properties[gNnode]
                    volume += sp.volume
                    occupancy_avg += sp.occupancy_avg
                    occupancy_std += sp.occupancy_std**2
                    center += sp.center
                N = float(len(gNnodes))
                if N > 0:
                    volume /= N
                    occupancy_avg /= N
                    occupancy_std = numpy.sqrt(occupancy_std)
                    center /= N
                    equivalence_name = np.equivalence_name  # redundant
                else:  # interstitial, just fill in bogus zero values
                    assert(np.label == SITELABEL['interstitial'])
                    center = None  # seems necessary so that the rec array is automatically built?!
                                   # (forces a |o4 list instead of 2D array)
                                   # (BUG: all fails if there's no interstitial!!!
                                   #   Need to construct rec array properly, define type, then fill array)
                    equivalence_name = 'interstitial'
                site_properties.append((np.label,volume,occupancy_avg,occupancy_std,center,equivalence_name,distance))
            # site_properties
            # (This is not the same site_properties that is used for
            # equivalence sites.)
            self.site_properties = numpy.rec.fromrecords(site_properties,
                                                         names='label,volume,occupancy_avg,occupancy_std,center,equivalence_name,distance')
            del site_properties
            self.site_properties.center[SITELABEL['interstitial']] = numpy.array([0.,0.,0.]) # hack for distance

            # build graph-indexed list of hopgraph-derived occupancies:
            # g_occupancy_avg, g_occupancy_std. shape = (order(2,combined_graph)) (cf g_label)
            self.g_occupancy_avg = numpy.zeros(self.g_labels.shape)
            self.g_occupancy_std = numpy.zeros(self.g_labels.shape)
            for ig,g in enumerate(self.graphs):
                nodes = numpy.array(g.graph.nodes())
                if not hasattr(g,'occupancy_avg'):
                    g.compute_site_occupancy()
                self.g_occupancy_avg[ig,self._graph2combined(ig,nodes)] = g.occupancy_avg[nodes]

            # compute distance of site from center of geometry of all sites (excl interstitial & bulk)
            centers = self.site_properties.center
            c0 = numpy.average(centers[SITELABEL['bulk']+1:])  # sites start after 'bulk'
            self.site_properties.distance = numpy.sqrt(numpy.sum((numpy.array(centers.tolist()) - c0)**2,axis=1))

            def _joined_dict(d, **kwargs):
                # needed below: build new dict from old dict and new key=value pairs
                _d = d.copy()
                _d.update(kwargs)
                return _d

            # add edges
            for ig,g in enumerate(self.graphs):
                # add graph number to properties to make edge unique
                # store the graph number in the NX 1.x edge dict as 'graph'
                ebunch = [(self._graph2combined(ig,u), self._graph2combined(ig,v),
                           _joined_dict(p, graph=ig))
                          for u,v,p in g.graph.edges(data=True)]
                # original pre NX 1.x implementation: p' <-- p + (ig,)
                # "(p'[3] == ig)"
                #ebunch = [(self._graph2combined(ig,u), self._graph2combined(ig,v),
                #          p+(ig,) )  for u,v,p in g.graph.edges]
                # Probably means that p[3] would always contain the graph number, p[:3]
                # would be other stuff ... k,N,fit ? -- OB 2010-06-28
                self.graph.add_edges_from(ebunch)
            # record common edges in edge_properties
            def edge_is_common(u,v):
                u0,v0 = self.g_labels[0,[u,v]]
                u1,v1 = self.g_labels[1,[u,v]]
                if not numpy.all([u0,v0,u1,v1]):  # get one or more None if node not common
                    return False
                return self.graphs[0].graph.has_edge(u0,v0) and self.graphs[1].graph.has_edge(u1,v1)
            for u,v,p in self.graph.edges(data=True):
                if self.edge_properties.has_key((u,v)):
                    self.edge_properties[(u,v)]['plist'].append(p['graph'])
                else:
                    self.edge_properties[(u,v)] = {'common': edge_is_common(u,v),
                                                   'plist':[p['graph']],}
        elif filename is not None:
            self.load(filename)
            return
        else:
            raise ValueError('Two HoppingGraphs or a filename is required.')

    def _graph2combined(self,igraph,x):
        """translate graph <igraph> node label(s) <x> to node label(s) of combined hybrid graph"""
        try:
            g2c = self._graph2combined_dict[igraph]
        except AttributeError:   # build cached dictionary first time
            self._graph2combined_dict = [None] * len(self.graphs)
            for ig in xrange(len(self.graphs)):
                 self._graph2combined_dict[ig] = dict(izip(self.g_labels[ig],self.node_properties.label))
            g2c = self._graph2combined_dict[igraph]
        if numpy.isscalar(x):
            return g2c[x]        # return scalar
        return numpy.array([g2c[label] for label in x])

    def _combined2graph(self,igraph,y):
        """translate  node label(s) <y> of combined hybrid graph to graph <igraph> node label(s)"""
        # or use self._graph2combined_dict[igraph] and invert... copy and paste for now:
        try:
            c2g = self._combined2graph_dict[igraph]
        except AttributeError:   # build cached dictionary first time
            self._combined2graph_dict = [None] * len(self.graphs)
            for ig in xrange(len(self.graphs)):
                 self._combined2graph_dict[ig] = dict(izip(self.node_properties.label,self.g_labels[ig]))
            c2g = self._combined2graph_dict[igraph]
        if numpy.isscalar(y):
            return c2g[y]        # return scalar
        return numpy.array([c2g[label] for label in y])

    def site_properties():
        doc = "site_properties of the combined graph, indexed by node label."
        def fget(self):
            return self.__site_properties
        def fset(self,x):
            self.__site_properties = x
        return locals()
    site_properties = property(**site_properties())

    def is_connected(self,igraph,n1,n2):
        """Return True if nodes n1 and n2 in graph igraph are connected."""
        graph = self.graphs[igraph]
        return graph.has_successor(n1,n2) or graph.has_predecessor(n2,n1)

    def load(self,filename=None):
        """Reinstantiate CombinedGraph from a pickled CombinedGraph (from save())."""
        filename = self.filename(filename,'pickle',set_default=True)
        fh = open(filename,'rb')
        try:
            h = cPickle.load(fh)
        finally:
            fh.close()
        # restore attributes from the temporary instance (works but not elegant...)
        self.graph = h.graph
        for attr in 'filtered_graph', 'graphs', 'g_labels', \
                'node_properties', 'edge_properties', \
                'g_occupancy_avg', 'g_occupancy_std':
            if hasattr(h,attr):
                self.__dict__[attr] = h.__dict__[attr]
        # managed attributes (evaluate for side effects)
        self.site_properties = h.site_properties
        del h

    def stats(self,igraph,data=None):
        """Statistics for the hopping graph.

        d = stats(igraph,[data=dict])

        Without the data argument, the method just returns some
        interesting values gathered from the graph igraph and the density. If
        a data dictionary is given, then the raw data are loaded into
        the dict and can be processed further by histogramming etc.

        :Arguments:
        igraph        number of the graph
        data          optional dictionary to hold raw data for
                      processing; modified by method

        :Returns:
        d             dictionary with expressive keys, holding the results
        """
        try:
            return self.graphs[igraph].stats(data=data)
        except IndexError:
            raise ValueError("igraph must be one of "+str(range(len(self.graphs))))


    def plot(self,igraph,filename=None,format='png',use_filtered_graph=True,label_sites=None,
             prog='neato', cmap=None, max_node_size=500, interactive=True, **drawargs):
        """Plot filtered graph using matplotlib.

        :Arguments:
        igraph         number of the graph (0 or 1)
        filename       file to write to
        format         any format that matplotlib allows and pdf
        use_filtered_graph
                       use a previously defined filtered graph (should be True)
        label_sites    {'all':False, 'common':True, 'none':False} switches that determine
                       which labels to add to the nodes
        prog           layout program, can be any of the graphviz programs
                       'dot','neato','twopi','circo','fdp','nop'
        cmap           matplotlib color map: nodes are colored by distance of the site from
                       the geometric center of all sites (excluding bulk)
        max_node_size  maximum node size (in point**2, q.v. matplotlib.scatter())
        interactive    True: display graph. False: only save to file (eg if no X11)
        **drawargs     additional keyword arguments to networkx.draw() (q.v.)
                       eg 'linewidths=(0.01,)' for vanishing outlines.
        """
        # Layout the combined graph and then plot each graph
        # separately but using the positions of the nodes of the
        # combined graph layout.
        import pylab
        from .utilities import matplotlib_interactive
        matplotlib_interactive(interactive)

        if label_sites is None:
            label_sites=dict(common=True,all=False,none=False)
        if cmap is None:
            cmap = dict(node=pylab.cm.jet,edge=None)

        graphviz_progs = ['dot','neato','twopi','circo','fdp','nop']
        graphnumbers = [0,1]
        if igraph not in graphnumbers:
            raise ValueError('Selected active graph <igraph> must be one of '+str(graphnumbers))
        if prog not in graphviz_progs:
            raise ValueError('Layout program <prog> must be one of '+str(graphviz_progs))
        graph = self.select_graph(use_filtered_graph)

        H = NX.MultiDiGraph(name=graph.name)  # need to clean up graph first --> H is tmp graph
        H.add_nodes_from([n for n in graph.nodes])
        # XXX: NX 1.x: original (u,v,{'graph':str(p[3])}) -- is the following correct, what was p[3] ??
        H.add_edges_from([(u,v,{'graph':str(p['graph'])}) for u,v,p in graph.edges(data=True)]) # clean up; str!

        pos = NX.pygraphviz_layout(H,prog=prog)  # combined graph layout

        def subgraph(ig,H=H):
            G = NX.DiGraph()
            G.add_nodes_from([n for n in H.nodes if self.node_properties[n].graphmask & 2**ig])
            G.add_edges_from([(u,v,p) for u,v,p in H.edges(data=True) if int(p['graph']) == ig])
            return G
        G = [subgraph(ig) for ig in graphnumbers]

        # labels
        if 'none' in label_sites and label_sites['none']:
            labels = dict(zip(graph.nodes(),[''] * graph.order()))
        else:
            labels = self.site_properties.equivalence_name.strip()         # sets dtype to |S12
            disjoint_sites = numpy.logical_not(self.node_properties.common)
            if 'all' in label_sites and label_sites['all']:
                labels[disjoint_sites] = self.node_properties.label[disjoint_sites]
            labels = dict(zip(graph.nodes(),labels[graph.nodes()]))   # only keep the ones in graph
        # colors
        vmin,vmax = self.site_properties.distance.min(),self.site_properties.distance.max()
        # shapes o = circle, d = diamond, s = square
        # only one per plot --> need to plot on top, not implemented yet
        ##common_nodes = [n for n in H.nodes if self.node_properties[n].common]
        ## #node_shape='o',nodelist=common_nodes,

        # for edge_color to work need networkx-0.37: https://networkx.lanl.gov/ticket/144
        # reorder graphnumbers so that active graph is last (and drawn on top)
        ordered_graphnumbers = graphnumbers[:]
        ordered_graphnumbers.remove(igraph)
        ordered_graphnumbers.append(igraph)
        # node size proportional to occupancy:
        # find normalization from all graphs (but exclude bulk)
        occ_max = []
        for g in self.graphs:
            try:
                occupancy = g.occupancy_avg
            except AttributeError:  # compute it
                g.compute_site_occupancy()
                occupancy = g.occupancy_avg
            realsites = numpy.ones(len(occupancy),dtype=numpy.bool)
            realsites[[SITELABEL['interstitial'],SITELABEL['bulk']]] = False
            occ_max.append(occupancy[realsites].max())  # masked out interstitial and bulk
        occupancy_max = numpy.max(occ_max)
        del occ_max

        # add something based on is_connected() to highlight sites NOT
        # connected to bulk, eg use different shape or heavier outline for
        # symbol
        #for ig,g in enumerate(self.graphs):
        #    g.has_bulkconnection = {}
        #    for n in g.nodes():
        #        g.has_bulkconnection[n] = self.is_connected(ig,n,SITELABEL['bulk'])

        def edgecolor(u,v):
            if self.edge_properties[(u,v)]['common']:
                return 'r'
            return 'k'
        def graph_parameters(ig,igraph=igraph):
            """Return dict of optional parameters for NX.draw.
            Extra function for interactive debugging."""
            if ig == igraph:
                alpha = 1
            else:
                alpha = 0.2
            node_color = self.site_properties.distance[G[ig].nodes()] # only keep the ones in graph
            gnodes = self._combined2graph(ig,G[ig].nodes())     # nodes in original graph
            occupancy = max_node_size * self.graphs[ig].occupancy_avg[gnodes]/occupancy_max
            edge_color = [edgecolor(u,v) for u,v in G[ig].edges]
            return dict(pos=pos,
                        node_color=node_color,
                        cmap=cmap['node'],
                        vmin=vmin,vmax=vmax,
                        node_size=occupancy,
                        edge_color=edge_color,
                        labels=labels,alpha=alpha,
                        )

        cax = pylab.axes(frameon=False)
        pylab.clf()
        for ig in ordered_graphnumbers:
            gp = graph_parameters(ig)
            gp.update(drawargs)
            NX.draw(G[ig],
                    # ax=cax,   ## needs to be commented out, otherwise empty frame
                    **gp)
        imgfile = self.filename(filename,format)
        pylab.savefig(imgfile)

    def export(self,igraph,filename=None,format='XGMML',imageformat=None,use_filtered_graph=True):
        """Layout the combined graph and highlight the chosen graph.

        h.export(igraph=0)

        :Arguments:
        graph        0 or 1, selects which graph is to be highlighted
        filename     name for the output files; appropriate suffixes are added
                     automatically
        format       XGMML or dot
        imageformat  graphics output format (png, jpg, ps, svg, ... see below)
        use_filtered_graph
                     By default, the filtered graph (see the filter() method) is
                     plotted. If set to False then the original HoppingGraph is
                     used instead.

        Common nodes are always highlighted in red and shown with the
        common label. Nodes and edges belonging to the selected graph
        are shown in black; the other graph is only shown in light
        gray.

        The graph is only written to an image file if an image format
        is supplied. See
        https://networkx.lanl.gov/reference/pygraphviz/pygraphviz.agraph.AGraph-class.html#draw
        for possible output formats but png, jpg, ps are safe bets.

        :Format:
        XGMML      http://www.cs.rpi.edu/~puninj/XGMML/draft-xgmml.html#Intro and
                   GML http://www.infosun.fim.uni-passau.de/Graphlet/GML/
        dot        See http://graphviz.org/doc/info/attrs.html for attributes.

        Note: On Mac OS X 10.3.9+fink the pygraphviz rendering is
        buggy and does not include node labels. Simply use the
        exported .dot file and use Mac OS X graphviz from
        http://www.pixelglow.com/graphviz/
        """
        # copy & paste and modified HoppingGraph.export() because I can't be bothered to write
        # a more generic version (or use site_properties properly?)
        import pygraphviz

        graphnumbers = range(len(self.graphs))
        if igraph not in graphnumbers:
            raise ValueError('igraph must be 0 or 1.')
        graph = self.select_graph(use_filtered_graph)
        if not format in self._export_formats:
            raise ValueError('Only the following formats are supported for export:'+
                             str(self._export_formats))
        activegraph = igraph            # index key for active graph
        graphnumbers.remove(igraph)
        inactivegraph = graphnumbers[0] # index key for inactive graph
        graphnumbers = range(len(self.graphs)) # need them all again

        G = pygraphviz.AGraph(directed=True)   # no multi edges??
        G.add_nodes_from(graph.nodes())
        G.add_edges_from([(e[0],e[1]) for e in graph.edges]) # n1 n2

        # attributes: http://www.graphviz.org/doc/info/attrs.html
        # colors: http://www.graphviz.org/doc/info/colors.html
        color = {'edge':     {'on' : 'gray20',
                              'off': 'gray75',
                              'common': 'red',
                              },
                 'nodefill': {'on' : 'gray60',
                              'off': 'gray90',
                              'common': 'black',
                              },
                 }
        # default to 'off' colours
        G.graph_attr['label'] = graph.name + ', k in 1/ns (N)'
        G.graph_attr['splines'] = 'true'
        G.node_attr['fontcolor'] = 'white'
        G.node_attr['shape'] = 'circle'
        G.node_attr['fillcolor'] = color['nodefill']['off']
        G.node_attr['style'] = 'filled'
        G.edge_attr['color'] = color['edge']['off']
        G.edge_attr['dir'] = 'forward'

        # colouring and labeling of edges
        # attach k_ji to edge (printed in 1/ns), number of transitions in ()
        xedgeattr = {}     # helper dict for XMMLG
        for e in G.edges():
            u,v = map(int,e)        # need nodes as integers
            xedgeattr[e] = {}            # hack for XMMLG
            xattr = xedgeattr[e]
            p0 = graph[u][v][0]          #  TODO: currently only chooses first edge 0 XXX
            k,N,fit,graph_number = p0['k'], p0['N'], p0['fit'], p0['graph'] # extract k from HoppingGraph
            if len(fit.parameters) == 3:         # get two rates if double exp
                a,k1,k2 = fit.parameters
                kns = 1000 * numpy.sort([k1,k2]) # rate in 1/ns
                k_dominant = kns[1]              # only display fastest rate
            else:
                k_dominant = 1000*k              # rate in 1/ns
            k_legend = "%(k_dominant).1f (%(N)d)" % locals()
            e.attr['headlabel'] = k_legend
            xattr['rate'] = k_dominant           # for XMMLG
            xattr['N'] = N                       # for XMMLG
            xattr['graphnumbers'] = [graph_number]  # for XMMLG, filter on graphnumber
            xattr['graph'] = graph_number
            xattr['common'] = False
            xattr['active'] = False
            if self.edge_properties[(u,v)]['common']:
                e.attr['color'] = color['edge']['common']
                xattr['graphnumbers'] = graphnumbers
                xattr['common'] = True
                xattr['active'] = True
            elif graph_number == igraph:
                e.attr['color'] = color['edge']['on']
                xattr['active'] = True

        # colouring and labeling of nodes
        xnodeattr = {}     # helper dict for XMMLG
        for k in G.nodes():
            node = G.get_node(k)
            label = int(node)
            common_label = self.node_properties[label].equivalence_name.strip()
            node_label = self.g_labels[activegraph,label]
            xnodeattr[label] = {}            # hack for XMMLG
            xattr = xnodeattr[label]
            xattr['common'] = False          # for XMMLG
            if common_label:
                node.attr['label'] = str(common_label)
                node.attr['fillcolor'] = color['nodefill']['common']
                xattr['graphnumbers'] = graphnumbers     # for XMMLG
                xattr['common'] = True                   # for XMMLG
                xattr['active'] = True                   # for XMMLG
                xattr['graph'] = activegraph
            elif node_label is not None:
                node.attr['label'] = str(node_label)
                node.attr['fillcolor'] = color['nodefill']['on']
                xattr['graphnumbers'] = [activegraph] # for XMMLG, filter on graphnumber (list with dot?)
                xattr['active'] = True                      # for XMMLG
                xattr['graph'] = activegraph
            else:
                node.attr['label'] = str(self.g_labels[inactivegraph,label])
                xattr['graphnumbers'] = [inactivegraph] # for XMMLG, filter on graphnumber (list with dot?)
                xattr['active'] = False                     # for XMMLG
                xattr['graph'] = inactivegraph

        if format == 'XGMML':
            filename = self.filename(filename,'.xgmml')
            xml = open(filename,'w')
            xml.write("""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n"""
                      """<!DOCTYPE graph PUBLIC "-//John Punin//DTD graph description//EN" "http://www.cs.rpi.edu/~puninj/XGMML/xgmml.dtd">\n"""
                      """<graph label="CombinedGraph" directed="1" id="CombinedGraph" """
                      """xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" """
                      """xmlns:dc="http://purl.org/dc/elements/1.1/" """
                      """xmlns:xlink="http://www.w3.org/1999/xlink" """
                      """xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" """
                      """xmlns="http://www.cs.rpi.edu/XGMML" """
                      """Vendor="hop.sitemap.CombinedGraph.export">\n"""
                      """\t<att name="documentVersion" value="1.0"/>\n""")
            xml.write("""\t<att name="activeGraph" type="integer" value="%d"/>\n""" % activegraph)
            for k in G.nodes():
                node = G.get_node(k)
                label = int(node)
                xattr = xnodeattr[label]
                xml.write("""\t<node id="%d" label="%s">\n""" % (label,node.attr['label']))
                xml.write("""\t\t<att type="list" name="graphs">\n""")
                for graphnumber in xattr['graphnumbers']:
                    xml.write("""\t\t\t<att type="integer" name="graphnumber" value="%d"/>\n""" \
                                  % graphnumber)
                xml.write("""\t\t</att>\n""")
                xml.write("""\t\t<att type="integer" name="graph" value="%(graph)d"/>\n""" % xattr)
                xml.write("""\t\t<att type="integer" name="common" value="%(common)d"/>\n""" % xattr)
                xml.write("""\t\t<att type="integer" name="active" value="%(active)d"/>\n""" % xattr)
#                xml.write("""\t\t<att type="" name="" value=""/>\n""")
                xml.write("""\t</node>\n""")
            for e in G.edges():
                u,v = map(int,e)
                xattr = xedgeattr[e]
                subst = {'u':u,'v':v,'label':e.attr['headlabel'],'id':"%d -> %d" % (u,v)}
                xml.write("""\t<edge source="%(u)d" target="%(v)d" label="%(label)s">\n""" \
                              % subst)
                xml.write("""\t\t<att type="string" name="EDGE_TYPE" value="DirectedEdge"/>\n""")
                xml.write("""\t\t<att type="string" name="interaction" value="DirectedEdge"/>\n""")
                xml.write("""\t\t<att type="integer" name="N" value="%(N)d"/>\n""" % xattr)
                xml.write("""\t\t<att type="real" name="rate" value="%(rate)f"/>\n""" % xattr)
                xml.write("""\t\t<att type="list" name="graphs">\n""")
                for graphnumber in xattr['graphnumbers']:
                    xml.write("""\t\t\t<att type="integer" name="graphnumber" value="%d"/>\n""" \
                                  % graphnumber)
                xml.write("""\t\t</att>\n""")
                xml.write("""\t\t<att type="integer" name="graph" value="%(graph)d"/>\n""" % xattr)
                xml.write("""\t\t<att type="integer" name="common" value="%(common)d"/>\n""" % xattr)
                xml.write("""\t\t<att type="integer" name="active" value="%(active)d"/>\n""" % xattr)
                xml.write("""\t</edge>\n""")
            xml.write("""</graph>\n""")
            xml.close()
        elif format == 'dot':
            dotname = self.filename(filename,'.dot')       # load into graphviz
            G.write(dotname)
            if imageformat:
                gfxname = filename+'.'+format
                G.layout(prog='dot')
                G.draw(path=gfxname,format=format)  # buggy on Mac OS X

    def plot_fits(self,**kwargs):
        raise NotImplementedError
    def export3D(self,**kwargs):
        raise NotImplementedError
    def tabulate_k(self,**kwargs):
        raise NotImplementedError
    def equivalent_sites_stats(self,graphnumber,elabels,equivalence=True):
        """Print statistics about one or a list of equivalence sites for the numbered graph.

        CombinedGraph.equivalent_sites_stats(graphnumber,elabels)

        :Arguments:
        graphnumber       index into CombinedGraph.graphs (typically, 0 or 1)
        elabels           single label or list of labels of equivalence sites
                          (without a '*' if the default identifier is used)
        equivalence       True: interprete elabels as equivalence labels
                          False: elabels are labels local to the graph (as used
                                 in the output of this method)
        """
        try:
            graph = self.graphs[graphnumber]
        except IndexError:
            raise ValueError('graphnumber must correspond to a graph in the CombinedGraph, one of '+
                             str(range(len(self.graphs))))
        graph.equivalent_sites_stats(elabels,equivalence=equivalence)

class fit_func(object):
    """Fit a function f to data (x,y) using the method of least squares.

    Attributes:

    parameters     list of parameters of the fit
    """
    # NOTE: this class must stay in the same file as HoppingGraph. Otherwise
    #       cPickle does not work on HoppingGraph (and one needs to implement
    #       custom pickling)
    def __init__(self,x,y):
        _x = numpy.asarray(x)
        _y = numpy.asarray(y)
        p0 = self.initial_values()
        fitfunc = self.f_factory()
        def errfunc(p,x,y):
            return  fitfunc(p,x) - y     # residuals
        p,msg = scipy.optimize.leastsq(errfunc,p0[:],args=(_x,_y))
        try:
            p[0]
            self.parameters = p
        except (TypeError,IndexError,):
            # TypeError for int p, IndexError for numpy scalar (new scipy)
            self.parameters = [p]
        self.message = msg

    def f_factory(self):
        """Stub for fit function factory, which returns the fit function.
        Override for derived classes.
        """
        def fitfunc(p,x):
            # return f(p,x); should be a numpy ufunc
            raise NotImplementedError("base class must be extended for each fit function")
        return fitfunc

    def initial_values(self):
        """List of initital guesses for all parameters p[]"""
        # return [1.0, 2.0, 0.5]
        raise NotImplementedError("base class must be extended for each fit function")

    def fit(self,x):
        """Applies the fit to all x values"""
        fitfunc = self.f_factory()
        return fitfunc(self.parameters,numpy.asarray(x))

class fitExp(fit_func):
    """y = f(x) = exp(-p[0]*x)"""
    def f_factory(self):
        def fitfunc(p,x):
            return numpy.exp(-p[0]*x)   # exp(-B*x)
        return fitfunc
    def initial_values(self):
        return [1.0]
    def __repr__(self):
        return "<fitExp "+str(self.parameters)+">"

class fitExp2(fit_func):
    """y = f(x) = p[0]*exp(-p[1]*x) + (1-p[0])*exp(-p[2]*x)"""
    def f_factory(self):
        def fitfunc(p,x):
            return p[0]*numpy.exp(-p[1]*x) + (1-p[0])*numpy.exp(-p[2]*x)
        return fitfunc
    def initial_values(self):
        return [0.5,0.1,1e-4]
    def __repr__(self):
        return "<fitExp2"+str(self.parameters)+">"

class fitlin(fit_func):
    """y = f(x) = p[0]*x + p[1]"""
    def f_factory(self):
        def fitfunc(p,x):
            return p[0]*x + p[1]
        return fitfunc
    def initial_values(self):
        return [1.0,0.0]
    def __repr__(self):
        return "<fitlin"+str(self.parameters)+">"


def survivalfunction(waitingtimes, block_w=200, block_t=1000):
    """Returns the survival function S(t), defined by a list of waiting times.

       survival([t0, t1, ..., tN]) --> S(t)

    S(t) is a function that gives the fractional number of particles that
    have not yet left the site after time t. It is 1 at t=0 and decays to 0.

    :Arguments:
       waitingtimes
          sequence of the waiting times from the simulations
      block_w
          reduce memory consumption by working on chunks of
          the waiting times of size <block_w>; reduce block_w if
          the code crashes with :Exception:`MemoryError`.
      block_t
          chunk input function arguments into blocks of size block_t

    :TODO: Make S(t) an interpolation function: massive speedup and fewer memory problems
    """
    S_doc = """Survival function S(t).  t can be a array.

    S(t) = <Theta(t0 - t)>_{t0} = <1-Theta(t-t0)>_{t0}
    """

    # NOTE: These big functions S_numpy or S_lowmem are *evaluated
    # every time* S(t) is called. Should be optimized eventually.

    # straightforward version
    def S_numpy(t):
        return numpy.sum([Unitstep(t0,t) for t0 in waitingtimes],axis=0)/len(waitingtimes)

    # work on chunks of waitingtimes of size block_w and chunks of times of size block_t
    def S_lowmem(t):
        _t = numpy.asarray(t)
        N_t = len(_t)
        func_values = numpy.zeros(N_t, dtype=float)    # final result S(t) = func_values[:]
        # need to block t, too, as it can be very large
        for j in numpy.arange(0, N_t, block_t):
            tmp_t = _t[j:j+block_t]
            values = numpy.zeros(len(tmp_t),dtype=float)
            for i in numpy.arange(0, len(waitingtimes), block_w):
                # the list comprehension can blow up in memory (so we block it)
                # TODO: rewrite this completely in a better way (cumsum?? Theta -> cumsum Delta)
                values += numpy.sum([Unitstep(t0,tmp_t) for t0 in waitingtimes[i:i+block_w]],axis=0)
            func_values[j:j+len(values)] = values
        func_values /= len(waitingtimes)
        return func_values

    # Perhaps only use S_lowmem ?? Would be more robust.
    if len(waitingtimes) < block_w:
        # print "survivalfunction(): using standard version (len(tau)=%d < block_w=%d)" % (len(waitingtimes),block_w)
        S = S_numpy
    else:
        # TODO: should really check if something like len(waitingtimes) * len(t) >
        # available memory but at this time t isn't known; I would have to put the
        # switching logic into the function S(t) itself The problem is really that for
        # long trajectories, len(t) can be very large (eg 100ns * 1000 frames/ns = 10^5)
        # when there are bound water molecules (with tau ~ t_sim); one way to solve it
        # would be to use a larger dt when evaluating S(t). Right now we are trying to be
        # somewhat robust by simply blocking both waitingtimes ('block_w') and times t
        # ('block_t')

        # print "survivalfunction(): using lowmem version (len(tau)=%d >= block_w=%d)" % (len(waitingtimes),block_w)
        S = S_lowmem
    S.__doc__ = S_doc
    return S

def Unitstep(x,x0):
    """Heaviside step function
                                       / 1    if x >= x0
    Unitstep(x,x0) == Theta(x - x0) = {  0.5  if x == x0
                                       \ 0    if x <  x0

    This is a numpy ufunc.

    :CAVEAT:    If both x and x0 are arrays of length > 1 then weird things are
                going to happen because of broadcasting.
                Using nD arrays can also lead to surprising results.

    :See also:  http://mathworld.wolfram.com/HeavisideStepFunction.html
    """
    _x = numpy.asarray(x)
    _x0 = numpy.asarray(x0)
    return 0.5*(1 + numpy.sign(_x - _x0))

