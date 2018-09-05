# setuptools installation of Hop
# Copyright (c) 2007-2013 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.
from __future__ import print_function

from setuptools import setup, find_packages
import versioneer

import sys
if sys.version_info[:2] < (2, 7) or sys.version_info[:1] >= (3,):
    print("HOP requires Python 2.7 (Python 3 not supported, sorry).  Python %d.%d detected" % \
        sys.version_info[:2])
    print("Please upgrade your version of python.")
    sys.exit(-1)

setup(name="Hop",
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description="Hop analyses solvent dynamics in molecular dynamics trajectories.",
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
      url="https://github.com/Becksteinlab/hop",
      keywords="science 'molecular dynamics' analysis hydration water",
      scripts = ['scripts/hop-generate-densities.py',
                 'scripts/hop-generate-hopgraph.py',
                 'scripts/hop-generate-hoptraj.py'],
      packages=find_packages(exclude=['tests','extras','doc/examples']),
      install_requires=['numpy>=1.0.3',
                        'scipy',
                        'networkx>=1.11',
                        'MDAnalysis>=0.18.0',
                        'GridDataFormats>=0.1.1',
                        ],
      extras_require={
          'plotting': ['matplotlib>=0.91.3',
                       'pygraphviz',         # only needed when plotting, not needed for graph building
                       ],
          'heatmap': ['rpy'],                # optional,used for heatmaps; or use rpy2
          },
)
