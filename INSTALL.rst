.. Hop Installation instructions

=========
 INSTALL
=========

This document gives an overview over additional requirements that
should be installed before Hop and Hop's installation process. The
latter should be painless because it only consists of pure python
scripts at the moment.

.. warning:: 
 
   This is legacy research code. It might not work at all. Use at your
   own risk. Feedback is very welcome through the `issue tracker`_.


.. _issue tracker: https://github.com/becksteinlab/hop/issues

.. Note:: Only Python 2.7 is currently supported.


.. _source-install:
	  
Source installation
-------------------

At the moment, only source installation is supported. Use pip_. Download the
tarball from https://github.com/Becksteinlab/hop/releases or do a
web-install:

.. code-block:: bash

   		pip install https://github.com/Becksteinlab/hop/archive/release-0.4.0.tar.gz

(This will also install all dependencies.)

.. _pip: https://pip.pypa.io


Conda
-----

At the moment, we do not have a conda package. However, it is easy to
set up a working environment for *hop* and the do the :ref:`install from
source <source-install>`_ described above.

.. code-block:: bash

		conda create -c conda-forge -n hop python=2.7 numpy scipy networkx MDAnalysis matplotlib pygraphviz GridDataFormats
		source activate hop

		pip install https://github.com/Becksteinlab/hop/archive/release-0.4.0.tar.gz
  

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

    
