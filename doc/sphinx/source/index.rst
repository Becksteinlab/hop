.. Hop documentation master file, created by
   sphinx-quickstart on Tue Jul 30 09:59:06 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Hop: analyzing solvent in molecular dynamics trajectories
=========================================================

:Release: |release|
:Date: |today|

This is a collection of python modules to analyze (primarily) water
behaviour in MD simulations. The idea is to find regions with a
density above a given threshold (hydration sites) and catalogue those
sites. Once this is done, one can analyze water movement in terms of
hops between those sites. The complicated solvation dynamics is thus
represented as a graph in which hydration sites are the nodes (or
vertices) and movements between sites are the edges.

Of course, it is also possible to look at the movement of other
particles such as ions or small molecules --- one simply selects a
different species.

The package is called *Hop* (no clever acronym, just quick to type,
and reflecting the fact that a "hopping analysis" is performed).

Hop is built with MDAnalysis_.

.. _MDAnalysis: https://www.mdanalysis.org


Contents
========

.. toctree::
   :numbered:
   :maxdepth: 1

   installation
   hop


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

