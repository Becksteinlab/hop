# $Id$
# setuptools installation of HOP
# Copyright (c) 2007-2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 2 (or higher, your choice)


from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup, find_packages

setup(name="HOP",
      version="0.1rc1",
      description="HOP analyses solvent dynamics in molecular dynamics trajectories",
      long_description="""\
      HOP performs a 'hopping analysis' of molecules in molecular
      dynamics (MD) trajectories. Typically, these molecules are water
      molecules. The movement of all waters is tracked as the move
      between hydration sites. Hydration sites are interpreted as
      vertices in a graph, while movement (or 'hops') between them are
      taken as edges. In this way hydration is characterized at a
      coarse grained level as a network of hops with rate constants
      and fluxes derived from the MD simulations.
      """,
      author="Oliver Beckstein",
      author_email="orbeckst@gmail.com",
      license="GPL2+",
      url="http://sbcb.bioch.ox.ac.uk/oliver/software/#HOP", # not set up yet
      keywords="science 'molecular dynamics' analysis hydration water",
      packages=find_packages(exclude=['tests','examples']),
      install_requires=['numpy>=1.0.3',
                        'networkx<0.99',       # **
                        'pygraphviz<0.99',     # **
                        'pylab>=0.91.3',
                        ],
      # also requires:  'MDAnalysis>0.5.1',
      # http://mdanalysis.googlecode.com/files/MDAnalysis-0.6.0-rc1.tar.gz
      # but this does not work automagically yet
      extras_require={'heatmap': ['rpy']},     # optional
      
)

#**) The code has not been updated yet to use the new NetworkX 1.0 API
#    hence only earlier versions (e.g. 0.36) are going to work reliably.

      
