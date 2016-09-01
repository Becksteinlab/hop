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
Warnings and exceptions
=======================

Custom warnings and exceptions.
"""


# Warnings
import warnings

class OverwriteWarning(Warning):
    """Warns that a variable or file will be overwritten or changed."""
    pass

class InconsistentDataWarning(Warning):
    """Warns that some input may not be consistent; in some cases it may
    actually be reasonable to supply such input and thus it does not raise an
    exception.

    If an exception is desired, use a warning filter, see
    http://docs.python.org/lib/warning-filter.html :

    >>> warnings.simplefilter('error',InconsistentDataWarning)
    """
    pass

class MissingDataWarning(Warning):
    """Warns that some data were not available; only a warning is raised
    because the code works around it or assumes default values (often 0).

    If an exception is desired, use a warning filter, see
    http://docs.python.org/lib/warning-filter.html :

    >>> warnings.simplefilter('error',MissingDataWarning)

    """
    pass

# Exceptions

class MissingDataError(Exception):
    """Signifies a critical error because required data is missing."""
    pass

class SelectionError(Exception):
    """Signifies an error with a atom selection."""
    pass

