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
Utility functions --- :mod:`hop.utilities`
==========================================

Random mix of convenience functions that don't fit anywhere else.

For messages I should probably use python's logger module but this is
working so far (even though it's pretty crappy).

"""

import sys
import os, errno
import cPickle
import warnings
import hop

def unlink_f(path):
    """Unlink path but do not complain if file does not exist."""
    try:
        os.unlink(path)
    except OSError, err:
        if err.errno != errno.ENOENT:
            raise

def mkdir_p(path):
    """Create a directory *path* with subdirs but do not complain if it exists.

    This is like GNU ``mkdir -p path``.
    """
    try:
        os.makedirs(path)
    except OSError, err:
        if err.errno != errno.EEXIST:
            raise


# unbound methods filename_function(), to be used in other
# classes; the plan is to make all classes that require them
# subclasses of hop.utilities.Saveable and bind them to this super
# class. (load() and save() are already tentatively XXXed out)

# used in many classes for filename handling (not all are Saveable yet)
#   filename = hop.utilities.filename_function
# (adds the _filename attribute to the class self!)
# NOTE: filename_function is not developed anymore and deprecated
#       Try to derive classes from Saveable.
def filename_function(self,filename=None,ext=None,set_default=False,use_my_ext=False):
    """Supply a file name for the object.

    fn = filename()             ---> <default_filename>
    fn = filename('name.ext')   ---> 'name'
    fn = filename(ext='pickle') ---> <default_filename>'.pickle'
    fn = filename('name.inp','pdf') --> 'name.pdf'
    fn = filename('foo.pdf',ext='png',use_my_ext=True) --> 'foo.pdf'

    The returned filename is stripped of the extension (use_my_ext=False) and
    if provided, another extension is appended. Chooses a default if no
    filename is given.  Raises a ValueError exception if no default file name
    is known.

    If set_default=True then the default filename is also set.

    use_my_ext=True lets the suffix of a provided filename take priority over a
    default ext(tension).
    """
    if filename is None:
        if not hasattr(self,'_filename'):
            self._filename = None        # add attribute to class
        if self._filename:
            filename = self._filename
        else:
            raise ValueError("A file name is required because no default file name was defined.")
        my_ext = None
    else:
        filename, my_ext = os.path.splitext(filename)
        if set_default:                  # replaces existing default file name
            self._filename = filename
    if my_ext and use_my_ext:
        ext = my_ext
    if ext is not None:
        if ext.startswith('.'):
            ext = ext[1:]  # strip a dot to avoid annoying mistakes
???LINES MISSING
???LINES MISSING
???LINES MISSING
???LINES MISSING
???LINES MISSING
???LINES MISSING
???LINES MISSING
???LINES MISSING
