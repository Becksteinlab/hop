# $Id$
# Hop --- a framework to analyze solvation dynamics from MD simulations
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
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

"""
Constants --- :mod:`hop.constants`
==================================

Constants that are being used throughout the hop module.

.. rubric:: Conversions:

The conversion factor f to a unit b' for a quantity X (whose numeric
value relative to the base unit b is stored in the program) is a
quantity with unit b'/b. In the dictionaries below only the numeric
value f(b->b') is stored::

  X/b' = f(b->b') * X/b

.. SeeAlso::

   :mod:`MDAnalysis.units`


Constants and functions
-----------------------

"""
from __future__ import absolute_import, division

import MDAnalysis.units
from MDAnalysis.units import (constants, water, conversion_factor, get_conversion_factor, convert)

#: Angstrom; conventional water van der Waals radius
r_water = 1.4

#: Avogradro's number in  mol**-1 from   http://physics.nist.gov/cgi-bin/cuu/Value?na
N_Avogadro = constants['N_Avogadro']


#: conversion_factor used by :func:`get_conversion_factor` (type, unit1, unit2):
#:
#: Note: any observable with a unit (i.e. one with an entry in
#:  Grid.unit) needs an entry in conversion_factor[] because
#: _check_set_unit() tests for type/value in this
#: dictionary. Derived classes can add to
#: conversion_factor, too, when they add new observable, eg threshold:
#: Another quantity with a unit (it is a density so it gets the
#: density unit factors)
conversion_factor['threshold'] = MDAnalysis.units.densityUnit_factor


#: pre-defined integer labels for sites --- DO NOT CHANGE THIS.
#: In particular, the order must be preserved and it is assumed that 'proper' sites
#: start after the bulk. The labels are also directly used as indices into arrays and
#: probably quite some code still implicitly assumes the definition below.
SITELABEL = dict(outlier = -1, interstitial = 0, bulk = 1)

