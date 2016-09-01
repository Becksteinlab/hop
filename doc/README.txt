===================
 Hop documentation
===================

:Author: Oliver Beckstein
:Date:   2010-10-28

The documentation for hop is a bit behind its development. If you are
a beta tester then you can probably save some time by just asking Oli
--- one can do many things with hop that are not immediately obvious.

A recent feature are simplified scripts (installed in the bin
directory and named ``hop-*.py``) that cover the basic three stages of
the hopping analysis.

1. generate densities (both solvent and bulk at the same time); note
   that at the moment the threshold is fixed but more options could be
   easily added to the script; it only depends on demand.

   required input
       RMS-fitted MD trajectory (anything that MDAnalysis_ can read)
       and a corresponding structure (i.e. topology file).
   output
       densities *water.pickle* and *bulk.pickle* (also dx files for
       visualization)

2. generate the hopping trajectory

   required input
        *water.pickle* and the same MD trajectory/structure as used
        for step 1
   output
        hopping trajectory *hoptraj.psf* and *hoptraj.dcd*

3. build the hopping graph, run basic analysis, export in various
   formats, and plot the rate constant data

   required input
        hopping trajectory *hoptraj.psf* + *hoptraj.dcd* and density
        *water.pickle*

   output
        hopping graph *hopgraph.pickle* and various other formats
        including *xgmml* for use in Cytoscape_; rate graphs are in
        the directory *survival_times*.


Any patches and improvements (including the docs) are very welcome.


.. _MDAnalysis: http://www.mdanalysis.org
.. _Cytoscape: http://www.cytoscape.org/
   
