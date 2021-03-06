========
 README
========

|docs| |zenodo|

**DEVELOPMENT VERSION of hop**: Please note that this is a beta
version of the package. It is still in heavy development. Feedback_ is
very welcome.


Hop: analyzing solvent in molecular dynamics trajectories
=========================================================

This is a collection of Python modules to analyze (primarily) water
behaviour in MD simulations. The idea is to find regions with a
density above a given threshold (hydration sites) and catalogue those
sites. Once this is done, one can analyze water movement in terms of
hops between those sites. The complicated solvation dynamics is thus
represented as a graph in which hydration sites are the nodes (or
vertices) and movements between sites are the edges.

Of course, it is also possible to look at the movement of other
particles such as ions or small molecules --- one simply selects a
different species.

The package is called 'Hop' (no clever acronym, just quick to type,
and reflecting the fact that a "hopping analysis" is performed).

Hop requires MDAnalysis_.

.. _MDAnalysis: https://www.mdanalysis.org


Installation
============

See the file `INSTALL.rst`_ for details.

The package is built on top of MDAnalysis_, which typically needs to be
installed beforehand.


.. _Install.rst:
   https://github.com/Becksteinlab/hop/blob/develop/INSTALL.rst


Documentation
=============

Please see the `online docs`_.

There is also some contents in the ``doc/`` directory, in particular
``doc/overview.txt``.

.. _online docs: https://hop.readthedocs.io


Bug reporting and Contributing
==============================

Please raise issues in the `issue tracker`_. This includes
problems or if you have suggestions how to improve the package or the
documentation. Pull requests for changes and fixes are also warmly
welcomed.

.. _issue tracker: https://github.com/Becksteinlab/hop/issues


Citing
======
|zenodo|

If you use Hop in published work please cite (for the time being) the
old abstract and the MDAnalysis paper (because Hop is built on top of
MDAnalysis):

* Ian Welland and Oliver Beckstein, (2015, June 23). hop: hop 0.3.4. Zenodo. http://doi.org/10.5281/zenodo.18864

* N Michaud-Agrawal, EJ Denning, TB Woolf, and O Beckstein. MDAnalysis: A Toolkit for the Analysis of Molecular
  Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,  doi:10.1002/jcc.21787

Thanks!

.. |docs| image:: https://readthedocs.org/projects/hop/badge/?version=latest
   :target: https://hop.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation
      
.. |zenodo| image:: https://zenodo.org/badge/13219/Becksteinlab/hop.svg
   :target: https://zenodo.org/badge/latestdoi/13219/Becksteinlab/hop
   :alt: Zenodo DOI

.. _Feedback: https://github.com/Becksteinlab/hop/issues
