# setuptools installation of Hop
# Copyright (c) 2007-2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup, find_packages

import sys
if sys.version_info[:2] < (2, 5):
    print "HOP requires Python 2.5 or better.  Python %d.%d detected" % \
        sys.version_info[:2]
    print "Please upgrade your version of python."
    sys.exit(-1)
if sys.version_info[:2] >= (2, 6):
    networkx_requirements = 'networkx>1.0'
else:
    # networkx 1.3 only works with 2.6+ so we fiddle the requirements
    networkx_requirements = 'networkx==1.2'

setup(name="Hop",
      version="0.3.3-devel",
      description="Hop analyses solvent dynamics in molecular dynamics trajectories",
      long_description="""\
Hop performs a 'hopping analysis' of molecules in molecular dynamics
(MD) trajectories. Typically, these molecules are water molecules. The
movement of all waters is tracked as they move between hydration
sites. Hydration sites are interpreted as vertices in a graph, while
movement (or 'hops') between them are taken as edges. In this way
hydration is characterized at a coarse grained level as a network of
hops with rate constants and fluxes derived from the MD simulations.\
""",
      author="Oliver Beckstein",
      author_email="orbeckst@gmail.com",
      license="GPLv3",
      url="http://sbcb.bioch.ox.ac.uk/oliver/software/#Hop",
      keywords="science 'molecular dynamics' analysis hydration water",
      scripts = ['scripts/hop-generate-densities.py',
                 'scripts/hop-generate-hopgraph.py',
                 'scripts/hop-generate-hoptraj.py'],
      packages=find_packages(exclude=['tests','extras','doc/examples']),
      package_data = {'vmd': ['*.tcl']},
      install_requires=['numpy>=1.0.3',
                        'scipy',
                        networkx_requirements,
                        'MDAnalysis>0.7.4',
                        'GridDataFormats>=0.1.1', # http://github.com/orbeckst/GridDataFormats
                        ],
      dependency_links = [
        "http://code.google.com/p/mdanalysis/downloads/list",
        "http://sbcb.bioch.ox.ac.uk/oliver/download/Python/",
        "http://github.com/orbeckst/GridDataFormats/downloads/",
        ],
      extras_require={
          'plotting': ['matplotlib>=0.91.3', # probably already installed
                       'pygraphviz',         # only needed when plotting, not needed for graph building
                       ],
          'heatmap': ['rpy'],                # optional,used for heatmaps; or use rpy2
          },
      zip_safe=True,          # vmdcontrol uses pkg_resources to find vmd tcl script
)
