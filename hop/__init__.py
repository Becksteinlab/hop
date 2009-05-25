# $Id$
# Copyright (c) 2007, 2008 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License.

__all__ = ['constants','sitemap','trajectory','graph','interactive',
           'utilities','external','analysis','siteanalysis','MCMC',
           'logger']

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

