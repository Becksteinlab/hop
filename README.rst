========
 README
========

|zenodo|


**DEVELOPMENT VERSION of hop**: Please note that this is a beta
version of the package. It is still in heavy development. Feedback is
very welcome (see below).


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

.. _MDAnalysis: http://www.mdanalysis.org/


Installation
============

See the file `INSTALL.rst`_ for details.

The package is built on top of MDAnalysis_, which typically needs to be
installed beforehand.


.. _Install.rst:
   https://github.com/Becksteinlab/hop/blob/master/INSTALL.rst


Documentation
=============

Please see the contents of the ``doc/`` directory, in particular
``doc/overview.txt``, and the Python doc strings.


Bug reporting
=============

Almost invariably, things will not work right away or it will be
unclear how to accomplish a certain task. In order to keep track of
feedback I ask you to use the Issue tracker at
http://github.com/Becksteinlab/hop/issues 

It helps tremendously to have everything in one place. Of course, feel
free to also additionally email me directly.

Thanks!


Citing
======
|zenodo|

If you use Hop in published work please cite (for the time being) the
old abstract and the MDAnalysis paper (because Hop is built on top of
MDAnalysis):

* Oliver Beckstein, Naveen Michaud-Agrawal and Thomas
  B. Woolf. Quantitative Analysis of Water Dynamics in and near
  Proteins. Biophysical Journal 96 (2009), 601a.
  doi:10.1016/j.bpj.2008.12.3147

* N Michaud-Agrawal, EJ Denning, TB Woolf, and O
  Beckstein. MDAnalysis: A Toolkit for the Analysis of Molecular
  Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
  doi:10.1002/jcc.21787

Thanks!



Contact
=======

Please do not hesitate to raise issues in the `issue tracker`_ or 
contact Oliver Beckstein <orbeckst@gmail.com> if problems occur 
or if you have suggestions how to improve the package or these instructions.

.. _issue tracker: http://github.com/Becksteinlab/hop/issues

.. |zenodo| image:: https://zenodo.org/badge/doi/10.5281/zenodo.18864.svg   
            :target: http://dx.doi.org/10.5281/zenodo.18864


