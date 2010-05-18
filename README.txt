.. $Id$

========
 README
========


Hop: analyzing solvent in molecular dynamics trajectories
=========================================================

This is a collection of python modules to analyze (primarily) water
behaviour in MD simulations. The idea is to find regions with a
density above a given threshold (hydration sites) and catalogue those
sites. Once this is done, one can analyze water movement in terms of
hops between those sites. The complicated solvation dynamics is thus
represented as a graph in which hydration sites are the nodes (or
vertices) and movements between sites are the edges.

The package is called 'hop' (no clever acronym, just quick to type).

Hop requires a bunch of other home-grown modules (notably my OpenDX
module) which are also included.


Installation
============

See INSTALL.txt for details.

The package makes heavy use of MDAnalysis
(http://mdanalysis.googlecode.com/).


Documentation
=============

Please see the contents of the doc/ directory, in particular
doc/overview.txt, and the python doc strings.


Contact
=======

Please do not hesitate to contact Oliver Beckstein
<orbeckst@gmail.com> if problems occur or if you have suggestions how
to improve the package or these instructions.


