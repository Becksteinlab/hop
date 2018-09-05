==============
 Installation
==============


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
web-install (choose the appropriate URL!):

.. code-block:: bash

   		pip install https://github.com/Becksteinlab/hop/archive/release-0.4.0.tar.gz

(This will also install all dependencies.)

.. _pip: https://pip.pypa.io


Conda
-----

At the moment, we do not have a conda_ package. However, it is easy to
set up a working environment for *hop* and then do the :ref:`install
from source <source-install>` described above.

.. code-block:: bash

		conda create -c conda-forge -n hop python=2.7 numpy scipy networkx MDAnalysis matplotlib pygraphviz GridDataFormats
		source activate hop
		pip install https://github.com/Becksteinlab/hop/archive/release-0.4.0.tar.gz

You can then run all the ``hop-*`` scripts in this environment.	(Exit
the environment as usual with ``source deactivate``.)

.. _conda: https://conda.io/docs/		
