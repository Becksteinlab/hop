.. Hop Installation instructions

=========
 INSTALL
=========

This document gives an overview over additional requirements that
should be installed before Hop and Hop's installation process. The
latter should be painless because it only consists of pure python
scripts at the moment.

Please do not hesitate to contact Oliver Beckstein
<orbeckst@gmail.com> if problems occur or if you have suggestions on
how to improve the package or these instructions.


Quick instructions
==================

If you are familiar with `pip`_ or `easy_install`_ from setuptools_ then try ::

  pip install Hop-0.3.3.tar.gz

or ::

  easy_install Hop-0.3.3.tar.gz

or ::

  cd Hop-0.3.3
  python setup.py install

Otherwise read the `Requirements`_ and `Installation`_ instructions below.

.. _pip: http://www.pip-installer.org/en/latest/
.. _`easy_install`: http://pythonhosted.org/setuptools/easy_install.html
.. _setuptools: http://pythonhosted.org/setuptools/

Requirements
============

If Hop is installed with `pip`_ or `easy_install`_ then some
dependencies are handled automatically; if not you will have to
install packages manually. This is always true for ``MDAnalysis`` at
the moment.

The ``GridDataFormats`` package can be obtained from my home page (see
below) but in most cases the installation should be able to
automatically download and install it.

Installing :mod:`matplotlib` (if plotting is desired) tends to be less
painful if it is done in advance, especially through your system's
package manager.


System requirements
-------------------

Tested with python 2.5, 2.6 (2.3 will not work anymore) on Linux and Mac OS X.


Required python packages
------------------------

.. Table:: Required packages for Hop.

   =============== ===================== ============================================================
   package         version               url
   =============== ===================== ============================================================
   MDAnalysis      >= 0.7.5              http://mdanalysis.googlecode.com
   numpy           >=1.0.3               http://numpy.scipy.org/
   NetworkX        >=1.0                 https://networkx.lanl.gov/
   GridDataFormats >= 0.1.1              https://github.com/orbeckst/GridDataFormats/archives/master
   =============== ===================== ============================================================

MDAnalysis_ requires additional modules; see the instructions at
its home page. Get the latest snapshot from
http://mdanalysis.googlecode.com/files/ or install it with ::

  pip install MDAnalysis

Note that installation of MDAnalysis_ can take a while because it
requires a python environment that is fully set up for scientific
computing. The good news is that once this is done then Hop should be
easy.

.. _MDAnalysis:: http://mdanalysis.googlecode.com


Optional packages
-----------------

The following packages should be installed to make best use of
Hop, especially for visualization.

.. Table:: Optional packages for Hop.

   =============== ===================== ============================================================
   package         version               url
   =============== ===================== ============================================================
   pygraphviz                            https://networkx.lanl.gov/pygraphviz/
   matplotlib       >=0.91.3             http://matplotlib.sourceforge.net/ 
   rpy (and R)                           http://rpy.sourceforge.net/
   biopython                             http://www.biopython.org
   =============== ===================== ============================================================

When these modules are missing the code will raise an :exc:`ImportError` when
attempting to use them; their absence is not always handled gracefuly
and it is suggested to install at least :mod:`matplotlib` (pylab). 

:program:`R` and :mod:`rpy` are only used for an experimental heatmap analysis
and can be safely ignored.

:mod:`pygraphviz` also requires :program:`graphviz`; see its own installation
instructions.



Optional external helper programs
---------------------------------

Some functions and classes can make use of external programs; none of
them are necessary for the core functionality of Hop but are listed
here for completeness.

VMD
     VMD's VolMap plugin can be used to generate densities. See
     `Visualization`_ for details on VMD.




Hints on obtaining packages
---------------------------

Many packages can be found through the local package manager (eg apt,
fink, yum, rpm). networkx is available at the above URL or with the
``easy_install`` command from python setuptools.

In Debian/Ubuntu::

   aptitude install python-setuptools pkg-config

   aptitude install graphviz graphviz-dev python-matplotlib
   easy_install networkx
   easy_install pygraphviz



Installation
============

Standard EasyInstall installation
---------------------------------

The Hop package can be installed with 'EasyInstall' which is part of
``setuptools`` (http://peak.telecommunity.com/DevCenter/setuptools).

If you have setuptools already installed you should be able to simply
install with ::

   easy_install Hop-0.3.3.tar.gz

You can install ``setuptools`` by following the instructions at
http://peak.telecommunity.com/DevCenter/EasyInstall#installation-instructions
or on Linux distributions such as Debian or Ubuntu by a
commandline similar to ::

   aptitude install python-setuptools

Note that ``easy_install`` tries to download additional packages via
the internet if it cannot find them installed already.


Installation in a custom location
---------------------------------

If you want to install Hop somewhere else but the system default
please have a look at
http://peak.telecommunity.com/DevCenter/EasyInstall#custom-installation-locations
which describes custom installations.

From python 2.6 onwards, an installation in a specific user directory
is supported. It is recommended to install in the following manner::

  cd Hop-0.3.3
  python setup.py install --user





Additional software
===================

The following software is not necessary to use the package itself but
has been found extremely useful by the author for using Hop or
analyzing data.


Interactive use and ``ipython``
-------------------------------

When Hop was developed, interactive use from a python command shell
turned out to be a very convenient application paradigm. ``ipython``
is very much recommended because of its ease to obtain interactive
help via ``?`` and ``??`` and to inspect objects via
TAB-completion. This is especially helpful because most of the
documentation is provided as python doc strings, both at the module
and at the class level.

For instance, to get an overview over interactive usage, load the
hop.interactive module and query the top level doc string::

 import hop.interactive
 hop.interactive ?


Visualization
-------------

VMD

  In order to visualize densities and water networks one can use `VMD
  <http://www.ks.uiuc.edu/Research/vmd/>`_, which reads natively psf
  and pdb files together with densities in OpenDX format.


Network analysis and layout
---------------------------

graphviz
  `graphviz <http://www.graphviz.org/>`_ is a versatile graph plotter,
  available for Linux and Mac OS X and integrated in most package
  systems. It is also required for ``pygraphviz``. It reads *dot*
  files.

Cytoscape
  `Cytoscape <http://www.cytoscape.org/>`_ is a very powerful network
  visualization platform; it reads *xgmml* files exported from Hop.

    
