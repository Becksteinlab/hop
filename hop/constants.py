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


Constants and functions
-----------------------

"""


#: Angstrom; conventional water van der Waals radius
r_water = 1.4

# all conversions: the conversion factor f to a unit b' for a
# quantity X (whose numeric value relative to the base unit b is
# stored in the program) is a quantity with unit b'/b. In the
# dictionaries below only the numeric value f(b->b') is stored.
#
#  X/b' = f(b->b') * X/b

#: conversion factors between the base unit (Angstrom) and other lengthUnits x:
#: L/x = L/A * lengthUnit_factor[x]
lengthUnit_factor = dict(Angstrom = 1.0, nm = 1.0/10, nanometer = 1.0/10)

# density conversion factor. Base unit is A**-3
# n/x = n/A**3 * densityUnit_factor[x]

# nm:
#   f = 1 A^-3/1 nm^-3 = 1/(10A)^-3 = 1/1000

# Molar:
#   factor = 1 A**-3 / (N_Avogadro * (10**-9 dm)**-3)

#: Avogradro's number in  mol**-1 from   http://physics.nist.gov/cgi-bin/cuu/Value?na
N_Avogadro = 6.02214179e+23

# relative to a density rho0 in g/cm^3:
#   M(H2O) = 18 g/mol   Molar mass of water
#
#   factor = 1/(1e-24 * N_Avogadro / M(H2O))
#     from rho/rho0 = n/(N_A * M**-1) / rho0  where [n] = 1/Volume, [rho] = mass/Volume

#: water density at T=298K, P=1atm
#:
#: .. Table: Water density from Jorgensen & Jenson, JCompChem 19 (1998), 1179
#:
#: ========= ===========
#:  model     g/cm^3
#: ========= ===========
#:   SPC     0.985(1)
#:   TIP3P   1.002(1)
#:   TIP4P   1.001(1)
#:   exp     0.997
#: ========= ===========
water = dict(exp=0.997, SPC=0.985, TIP3P=1.002, TIP4P=1.001)  # in g cm**-3
water['MolarMass'] = 18.016                                   # in g mol**-1

#: for water densities, this is the volume per water molecule in A**3
densityUnit_factor = dict(
    Angstrom=1/1.0,
    nm=1/1e-3, nanometer=1/1e-3,
    Molar = 1/(1e-27*N_Avogadro),
    SPC = 1/(1e-24*N_Avogadro*water['SPC']    / water['MolarMass']),
    TIP3P = 1/(1e-24*N_Avogadro*water['TIP3P']/ water['MolarMass']),
    TIP4P = 1/(1e-24*N_Avogadro*water['TIP4P']/ water['MolarMass']),
    water = 1/(1e-24*N_Avogadro*water['exp']  / water['MolarMass']),
    )


#: basic time unit is ps;
#: 1 AKMA time unit = 4.888821E-14 sec
#:    http://brooks.scripps.edu/charmm_docs/c28docs/c28b2/html/usage.html#AKMA
#:  1200ps/ps * f = 1.2 ns/ns  ==> f = 1/1000
timeUnit_factor = dict(
    ps=1/1.0,
    ns=1/1e3,
    second=1/1e12,
    AKMA=1/4.888821e-2,
    )

#: conversion_factor used by :func:`get_conversion_factor` (type, unit1, unit2):
#:
#: Note: any observable with a unit (i.e. one with an entry in
#:  Grid.unit) needs an entry in conversion_factor[] because
#: _check_set_unit() tests for type/value in this
#: dictionary. Derived classes can add to
#: conversion_factor, too, when they add new observable, eg threshold:
#: Another quantity with a unit (it is a density so it gets the
#: density unit factors)
conversion_factor = dict(length=lengthUnit_factor,
                         density=densityUnit_factor,
                         threshold = densityUnit_factor,
                         time=timeUnit_factor,
                         )


#: pre-defined integer labels for sites --- DO NOT CHANGE THIS.
#: In particular, the order must be preserved and it is assumed that 'proper' sites
#: start after the bulk. The labels are also directly used as indices into arrays and
#: probably quite some code still implicitly assumes the definition below.
SITELABEL = dict(outlier = -1, interstitial = 0, bulk = 1)

def get_conversion_factor(unit_type,u1,u2):
    """generate the conversion factor u1 -> u2 by using the base
    unit as an intermediate

    conversion_factor = get_conversion_factor('density','SPC','Molar')

    f[u1 -> u2] = factor[u2]/factor[u1]

    Conversion of X (in u1) to X' (in u2):

        X' = conversion_factor * X
    """
    # x is in u1: from u1 to b:  x'  = x  / factor[u1]
    #             from b  to u2: x'' = x' * factor[u2]
    # so f[u1,u2] = factor[u2]/factor[u1]
    return conversion_factor[unit_type][u2] / conversion_factor[unit_type][u1]
