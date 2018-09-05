# Hop --- a framework to analyze solvation dynamics from MD simulations
# Copyright (c) 2009-2014 Oliver Beckstein <orbeckst@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""\
Hop --- a framework to analyze solvation dynamics from MD simulations
=====================================================================

This is a collection of python modules to analyze (primarily) water
behaviour in molecular dynamics (MD) simulations. The idea is to first
find regions with a density above a given threshold (hydration sites)
and catalogue those sites. Once this is done, one can analyze water
movement in terms of hops between those sites. The complicated
solvation dynamics is thus represented as a graph in which hydration
sites are the nodes (or vertices) and movements between sites are the
edges.

Please see the python doc strings and the documentation in the doc/
directory. The :mod:`hop.interactive` module contains convenience wrapper
functions primarily for interactive use in :program:`ipython` but it also has
extensive documentation and the functions act as examples for how to
use the module.

>>> import hop.interactive
>>> help(hop.interactive)

"""

__all__ = ['constants','sitemap','trajectory','graph','interactive',
           'utilities','analysis','siteanalysis','MCMC',
           'qhull']

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
