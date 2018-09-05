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
Quickstart: using the hop package --- :mod:`hop.interactive`
============================================================

A typical session starts with a trajectory (which should have been
RMS-fitted to a reference structure). Any topology and trajectory file
suitable for MDAnalysis_ can be used such as PSF+DCD, PDB+XTC or a
single PDB. In the following Charmm/NAMD psf and dcd files are used as
examples.

We will use the high-level wrapper functions in :mod:`hop.interactive`:

>>> from hop.interactive import *

.. _MDAnalysis:: http://mdanalysis.googlecode.com

Hydration sites
---------------

Hydration sites are sites of water density higher than the bulk
density but one special site is the bulk. The hydration sites and the
bulk site are computed in two separate steps.


High density sites
..................

First build the density of the water oxygens.

>>> density = make_density(psf,dcd,filename,delta=1.0)

If you have VMD with the VolMap plugin installed *and* your trajectory
fits into your computer's RAM you can also choose the VMD backend to
compute the density (which can be marginally faster):

>>> density = make_density(psf,dcd,filename,delta=1.0,backend='VMD')


The density is also saved as a pickled python object so that one can
easily reload it. The density is also exported as a dx file for
visualization (e.g. use :func:`hop.interactive.visualize_density`,
which calls :program:`VMD`).

From the density one creates the 'site map' for a given threshold (by
default this is a multiple of the water bulk density):

>>> density.map_sites(threshold=2.72)

Experiment with the threshold; :class:`hop.analysis.DensityAnalysis` can help
to systematically explore the parameter space, and it is also helpful
to load the density into a visualization software such as VMD and
interactively explore contour levels. Values between 1.65 and 3 have
given decent results in the past but this is system-dependent.)


Bulk site
.........

For a full analysis of hopping events one also needs to define a bulk
site. This is currently accomplished by calculating a second 'bulk'
density (all water not within 3.5 A of the protein) and manually
inserting the bulk site into the site map for the first density.

>>> density_bulk = make_density(psf,dcd,'bulk',delta=1.0,
...         atomselection='name OH2',
...         soluteselection='protein and not name H*',
...         cutoff=3.5
...         )

Using VMD's VolMap can be potentially be faster --- try it if the
default seems too slow to you:

>>> density_bulk = make_density(psf,dcd,'bulk',delta=1.0,
...         atomselection='name OH2 and not within 3.5 of (protein and name not hydrogen)',
...         backend='VMD',load_new=False)

The bulk density should be a big, well defined volume so we choose a
fairly low threshold:

>>> density_bulk.map_sites(0.6)

Add the biggest bulk site:

>>> density.site_insert_bulk(density_bulk)
>>> density.save()
>>> del density_bulk

.. Note:: Behind the scenes, the bulk is simply prepended to the list
   of all sites (``density.sites``,
   :attr:`hop.sitemap.Density.sites`), found so far. By convention the
   site at position 1 in the list of all sites is treated specially in
   many parts of hop (it has the so-called sitelabel "1", which is
   simply the position in the list of sites) and hence you might
   encounter unexpected behaviour later if you do not insert a bulk
   site.

Statistics about the sites can be produced with

>>> analyze_density(density,figname)

The results figures will be named `<figname>.pdf`.


Remapping for comparing site maps
.................................

This section is only relevant if you plan on comparing site maps. Then
you *must* compare the density to your reference density *now* before
proceeding. You will

   1) remap this density to be defined on the same grid as the reference
      density (for this to work, this density must have been generated from
      a trajectory that has been RMS-fitted to the same reference structure
      as; see hop.trajectory.rms_fit_trj() and
      hop.trajectory.fasta2select())

      >>> ref_density = hop.sitemap.Density(filename='my_reference_density')
      >>> remapped_density = hop.sitemap.remap_density(density,ref_density)

   2) find the equivalence sites in the two densities and add those sites
      to **both** densities:

      >>> remapped_density.find_equivalence_sites_with(ref_density,verbosity=3)
      >>> remapped_density.save(<filename>)
      >>> ref_density.save()

      (You must also recalculate the reference densities hopping
      trajectory (see below) because some sites may have been merged into
      'equivalence sites'. See docs for
      hop.sitemap.find_equivalence_sites_with() and
      hop.graph.CombinedGraph()).

      From now on, work with the remapped density:
      >>> density = remapped_density


Hopping trajectory
------------------

Next we translate the dcd into a 'hopping trajectory' (saved in dcd
format) in which coordinates for a given water oxygen are replaced by
the site it visits at each time step.

>>> hops = make_hoppingtraj(density,'hop_water+bulk')

All further analysis should use this hopping trajectory (from disk) as
it is computationally much cheaper to read the trajectory than to
re-translate the coordinate trajectory (which is done behind the
scences if the hopping trajectory is not available).


Hopping graph
-------------

The final step is to map out the graph of transitions between sites
(using the hopping trajectory):

>>> tn = build_hoppinggraph(hops,density)

tn.hopgraph holds this graph (tn.graph just contains all jumps
including the interstitial and off-sites). The edges of hopgraph are
the rate constants k_ji (in 1/ps) for hops i --> j. They are computed
from an exponential fit to the site survival function S_ji(t) for
particles waiting to hop from i to j.

The density is provided to attach data to the nodes of the
hopgraph. It is required for visualization and analysis (although not
strictly necessary for the hopgraph itself).

Further analysis uses tn.hopgraph:

>>> h = tn.hopgraph           # main result is the 'hopgraph'
>>> h.save('hopgraph')        # save the hopping graph (necessary for cg part)
>>> h.filter(exclude={'outliers':True, 'Nmin':2, 'unconnected':True})
>>> h.show_rates()            # show all calculated rate constants (filtered graph)
>>> h.plot_fits(xrange(301))  # plot rate constant fits for t=0ps to 300ps
>>> h.plot_fits()
>>> h.export('water')         # write dot file to visualize (filtered) graph
>>> h.plot_site_occupancy('siteoccupancy')  # plot site occupancy from graph ---is NOT working
>>> h.plot_residency_times('residencytimes')# residency times --- is NOT working

To compare the water network based on density with another hop graph
(based on ref_density), construct the CombinedGraph:

>>> h_ref = hop.graph.HoppingGraph(filename=<filename>) --- basically repeat steps from
###                                                     --- ref_density only with differ labels
>>> cg = hop.graph.CombinedGraph(g0=h,g1=h_ref)
>>> cg.plot(0,'cg_h',linewidths=(0.01,))
>>> cg.plot(1,'cg_h_ref',linewidths=(0.01,))


.. rubric:: TODO

Currently un(der)-documented:
* Remapping densities to a reference density (see hop.sitemap.remap_density).
* Comparing densities and finding equivalence sites (see
  hop.sitemap.find_common_sites() and Density.find_equivalence_sites_with()).
* Comparing hopgraphs across different simulations: requires equivalence sites in
  both densities; then build the hop.graph.CombinedGraph().

Functions
---------

"""
from __future__ import absolute_import

import os
import logging

import numpy
import MDAnalysis

from . import sitemap
from . import density
from . import trajectory
from . import graph
from . import constants
from .density import DensityCreator

logger = logging.getLogger("hop.interactive")

def generate_densities(*args, **kwargs):
    """Analyze the trajectory and generate solvent and bulk density.

    generate_densities(topol, traj, atomselection='name OW') --> densities

    This function can take a long time because it has to read the whole
    trajectory. Progress is printed to the screen. It saves results to pickle
    files. These files are :class:`hop.sitemap.Density` objects and can be used
    to instantiate such a density object.

    :Arguments:
      filename
         name of the solvent density with bulk site
      bulkname
         bulk density
      density_unit
         unit of measurement for densities and thresholds
         (Molar, nm^{-3}, Angstrom^{-3}, water, SPC, TIP3P, TIP4P)
      solvent_threshold : exp(1) = 2.7182818284590451
         hydration sites when density > this threshold
      bulk_threshold : exp(-0.5) = 0.60653065971263342
         bulk site are regions with density > this threshold
         (and water farther away from the protein heavy atoms than *cutoff*)
      delta : 1.0
         cubic grid size in Angstrom
      cutoff
         bulk-water is assumed to start at this distance from the
         soluteselection
      soluteselection : "protein and not name H*"
         how to select the solute (for bulk density)

    :Returns: a dict containing :class:`hop.sitemap.Density` instances for the
              the "solvent" and the "bulk" density; the "solvent" has the bulk
              site (largest site in "bulk") inserted as site 1.

    .. Note::

       The "solvent" density is going to be used throughout the rest of the
       protocol. Should you ever remap the sites (i.e. run
       :meth:`~hop.sitemap.Density.map_sites` with a different threshold) then
       **you must insert the bulk site again** (because the bulk site is
       removed for technical reasons whenever the sites change); use the saved
       bulk site and the :meth:`hop.sitemap.Density.site_insert_bulk` method.

    .. SeeAlso::

       Keyword arguments are passed on to :class:`hop.density.DensityCreator`
       where all possible keywords are documented; the site mapping is done
       with :meth:`hop.sitemap.Density.map_sites`.
    """
    filename = kwargs.pop('filename', 'water')  # solvent pickle file
    bulkname = kwargs.pop('bulkname', 'bulk')   # bulk solvent pickle file
    solvent_threshold = kwargs.pop('solvent_threshold', numpy.e)
    bulk_threshold = kwargs.pop('bulk_threshold', numpy.exp(-0.5))
    density_unit = kwargs.pop('density_unit', "water")
    kwargs['mode'] = "all"
    DC = DensityCreator(*args, **kwargs)
    densities = DC.create()
    # save the precious files right away
    densities['bulk'].save(bulkname)
    densities['solvent'].save(filename)
    # modifies densities['solvent'] (and also returns it)
    density = DC.DensityWithBulk(density_unit=density_unit,
                                 solvent_threshold=solvent_threshold,
                                 bulk_threshold=bulk_threshold)
    # save again, but now with mapped sites and bulk site include in "solvent"
    densities['bulk'].save(bulkname)
    density.save(filename)  # This is solvent but with bulk site added.
    # finally write dx files for visualization
    densities['bulk'].export()
    density.export()
    logger.info("generate_densities(): saved densities and exported dx files")
    return densities

def make_density(psf,dcd,filename,delta=1.0,atomselection='name OH2',
                 **kwargs):
    """Build the density by histogramming all the water oxygens in a dcd.

    density = make_density(psf,dcd,filename,delta=1.0)

    The function builds the density object, writes it to disk, and
    also exports it as a dx file for visualization (use
    vizualize_density(density)).

    :Arguments:

      *psf
         topology
      *dcd*
         trajectory (should be RMS fitted to a reference frame)
      *filename*
         default filename for the density
      *delta*
         grid spacing Angstrom

      *kwargs*:
            *padding*
               increase box dimensions for 3D histogramming by padding
            *soluteselection*
            *cutoff*
                for bulk density: setting both `soluteselection='protein and not name H*'`
                and `cutoff=3.5` A selects *'<atomsel> NOT WITHIN <cutoff> OF <solutesel>'*

    :Returns: *density*, :class:`hop.sitemap.Density` object; the density is
              converted to a fraction of the density of bulk TIP3P water

    """
    density = density.density_from_trajectory(
        psf,dcd,delta=delta,atomselection=atomselection,
        verbosity=3,**kwargs)
    density.convert_density('TIP3P')
    density.save(filename)
    density.export()

    logger.info("make_density() completed. Run map_sites(density,...) next")
    return density

def make_xstal_density(pdb,filename,**kwargs):
    """Generate a density from the crystalwaters in a PDB.

    For arguments see :func:`make_density`.

    (These water are typically named HOH.)

    .. SeeAlso:: Water molecules are counted as point-like
       particles. One can also use
       :class:`hop.sitemap.BfactorDensityCreator` to broaden water
       molecules according to their B-factor.

    """
    kwargs.setdefault('atomselection', "resname HOH")
    return make_density(pdb,filename,**kwargs)

def analyze_density(density,figure='sitestats'):
    """Site statistics based on the density alone.

    Plots site volumes, average densities and occupancy, and writes it to the
    pdf file <figure>.pdf
    """
    import pylab
    #import os

    # convert density to the chosen density unit (typically, relative to bulk water)
    factor = constants.get_conversion_factor('density',
                                             density.unit['length'] + "^{-3}",
                                             density.unit['density'])

    x,N,DN = density.site_occupancy(include='sites')
    x,V = density.site_volume(include='sites')

    pylab.clf()
    pylab.title('Site statistics for threshold %.2f (%s)' %
                (density.P['threshold'],density.unit['threshold']))
    pylab.plot(x,V/10,'ro')
    pylab.plot(x,factor * N/V,'r.')
    pylab.plot(x,N,'b.')
    pylab.legend(('V/10 %s**3' % density.unit['length'],
                  'density %s' % density.unit['density'],
                  'occupancy'))
    pylab.errorbar(x,factor*N/V,yerr=factor*DN/V,fmt=None,ecolor='r',capsize=0)
    pylab.errorbar(x,N,yerr=DN,fmt=None,ecolor='b',capsize=0)
    pylab.xlabel('site label')

    #eps = figure+'.eps'
    pdf = figure+'.pdf'
    pylab.savefig(pdf)
    #os.system('epstopdf --outfile=%(pdf)s %(eps)s' % locals())
    #os.remove(eps)

def visualize_density(density):
    """Visualize the trajectory with the density in VMD.

       visualize_density(density)

    :Arguments:
       *density*
          hop.sitemap.Density object

    """
    dx = density.filename() + '.dx'
    os.system('vmd '+density.metadata['psf']+' '+density.metadata['dcd']+' -m '+dx)

def make_hoppingtraj(density,filename,**hopargs):
    """Create the hopping trajectory from a density with a site map.

       hops = make_hoptraj(density,filename)

    :Arguments:
      *density*
         density object with a site map
      *filename*
         prefix for the hop trajectory files (psf and dcd)
      *hopargs*
         keyword args to add to :class:`~hop.trajectory.HoppingTrajectory` such
         as `fixtrajectory = {'delta':10.22741474887299}`

    This function relies on the density's metadata. In particular it
    uses `density.metadata['psf']` and `metadata['dcd']` to find its input
    data and `metadata['atomselection']` to define the atoms to track.
    """
    try:
        if len(density.sites) < 2:
            raise ValueError
    except AttributeError,ValueError:
        errmsg = 'The density misses a site map or has only one site.'
        logger.fatal(errmsg)
        raise ValueError(errmsg)

    u = MDAnalysis.Universe(density.metadata['psf'],density.metadata['dcd'])
    group = u.select_atoms(density.metadata['atomselection'])
    hops = trajectory.HoppingTrajectory(u.trajectory,group,density,**hopargs)
    hops.write(filename)
    return hops

def build_hoppinggraph(hoppingtrajectory,density):
    """Compute the graph of all site hops and calculate the rate constants.

    tgraph = build_hoppinggraph(hops,density)

    :Arguments:
    hops
           hop.trajectory.HoppingTrajectory object
    density
           hop.sitemap.Density object

    :Returns:  tgraph, a :class:`hop.graph.TransportNetwork` object
    """
    tgraph = graph.TransportNetwork(hoppingtrajectory,density)
    tgraph.compute_residency_times()
    return tgraph

def build_hoppinggraph_fromfiles(hoppingtrajectory_filename,density_filename):
    """Compute the TransportNetwork including HoppingGraph from files.

    tn = build_hoppinggraph_fromfiles('hoptraj','water_density')

    Input:
    hoppingtrajectory_filename   filename for HoppingTrajectory psf and dcd
    density_filename             filename for pickled Density

    Output:
    tn                           hop.graph.TransportNetwork object (qv)
    """
    hoppingtrajectory = trajectory.HoppingTrajectory(filename=hoppingtrajectory_filename)
    density = sitemap.Density(filename=density_filename)
    return build_hoppinggraph(hoppingtrajectory,density)

def hopgraph_basic_analysis(h, density, filename):
    """Do some simple analysis tasks on the hopgraph.

    hopgraph_basic_analysis(h, density, filename)

    :Arguments:
       h
          hopgraph, a :class:`hop.graph.HoppingGraph`
       density
          density, a :class:`hop.sitemap.Density`
       filename
          default filename for generated files; all files and new
          directories are written in the directory pointed to by the
          path component
    """
    analysisdir = os.path.dirname(filename)
    logger.warn("Setting analysisdir = %(analysisdir)r", vars())

    ratesfile = os.path.join(analysisdir, 'rates.txt')
    h.show_rates(filename=ratesfile)
    logger.info("Wrote all rates to %(ratesfile)r.", vars())

    h.filter(exclude={'outliers':True, 'bulk':True})
    h.export(filename, format='XGMML')
    logger.info("Exported hopgraph as %(filename)s.xgmml", vars())

    logger.info("Generating 3D graph %(filename)s.psf/pdb", vars())
    logger.info("Note: bulk site omitted for clarity.")
    h.export3D(density)

    survival_time_dir = os.path.join(analysisdir, 'survival_times')
    logger.info("Generating survival times plots in %(survival_time_dir)r", vars())
    logger.info("This takes a while. Note: transitions to/from bulk are excluded, Nmin=5.")

    h.filter(exclude={'outliers':True, 'bulk':True, 'Nmin':5})
    h.plot_fits(directory=survival_time_dir, ncol=2, nrow=3)
