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

.. warning:: This is legacy software that is provided "AS IS". In
	     particular, there are currently *no tests* and it is not
	     guaranteed to work or produce correct results. Help and
	     contributions are welcome!


License
-------

*hop* is released under the `GNU General Public License, v3`_ (because
it links to MDAnalysis_, which is GPL licensed).	  

.. _GNU General Public License, v3: 
   https://www.gnu.org/licenses/gpl-3.0.en.html


Documentation
-------------

The primary documentation consists of the `online docs`_ (which you
are reading).

There is also some content in the ``doc/`` directory, in particular
``doc/overview.txt``.

.. _online docs: https://hop.readthedocs.io


Bug reporting
-------------

Almost invariably, things will not work right away or it will be
unclear how to accomplish a certain task. In order to keep track of
feedback I ask you to use the Issue tracker at https://github.com/Becksteinlab/hop/issues


Citing
------

If you use Hop in published work please cite (for the time being) the abstract
[Hop2009]_ and the MDAnalysis paper (because Hop is built on top of MDAnalysis)
[MDAnalysis2011]_. Thanks!

.. [Hop2009] Oliver Beckstein, Naveen Michaud-Agrawal and Thomas B.
   Woolf. Quantitative Analysis of Water Dynamics in and near
   Proteins. Biophysical Journal 96 (2009), 601a.
   doi:10.1016/j.bpj.2008.12.3147

.. [MDAnalysis2011] N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein. 
   MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics
   Simulations. J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787


Contents
--------

.. toctree::
   :numbered:
   :maxdepth: 1

   installation
   hop


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

