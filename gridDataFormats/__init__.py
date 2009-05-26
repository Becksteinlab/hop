# $Id$
# gridDataFormats --- python modules to read and write gridded data
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Lesser Public License, version 3 or later.

"""gridDataFormat module

Copyright (c) 2007, 2008 Oliver Beckstein <orbeckst@gmail.com>
Released under the GNU Public License.

This module contains classes that allow importing and exporting of
simple gridded data, A grid is an N-dimensional array that represents
a discrete mesh over a region of space. The array axes are taken to be
parallel to the cartesian axes of this space. Together with this array
we also store the edges, which are are (essentially) the cartesian
coordinates of the intersections of the grid (mesh) lines on the
axes. In this way the grid is anchored in space.

The Grid class acts as a universal constructor for specific formats.

 g = Grid(**kwargs)           # construct
 g.export(filename, format)   # export to the desire format

Some formats can also be read:

 g = Grid()                   # make an empty Grid
 g.load(filename)             # populate with data from filename

See the doc string for Grid for details on **kwargs.

Formats:

   OpenDX        IBM's Data Explorer, http://www.opendx.org/
   gOpenMol      http://www.csc.fi/gopenmol/
   pickle        python pickle file
"""

__all__ =  ['core','OpenDX','gOpenMol']

import warnings

class gridDataWarning(Warning):
    """Warns of a problem specific to the gridData module."""
    pass

from core import Grid
