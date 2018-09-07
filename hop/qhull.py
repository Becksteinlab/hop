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
Using qhull to define regions for hopping analysis --- :mod:`hop.qhull`
=======================================================================

Interface to some functions of the `qhull`_ (or rather the `qconvex`_)
program. `qhull`_ must be installed separately (see links).

The main functionality is to define a region in space within the convex hull of
a protein. The hull is typically defined by a selection of atoms and written as
a "density" file for use in :mod:`hop`.

.. _qhull: http://www.qhull.org/html/index.htm
.. _qconvex: http://www.qhull.org/html/qconvex.htm

Example
-------

In this example the convex hull of the C-alpha atoms is
computed. Initially, points must be extracted from the structure to a file::

  hop.qhull.points_from_selection(psf='protein.psf', pdb='protein.pdb', filename='ca_100%.dat')

and saved to file ``ca_100%.dat``.

This is usually too large and also entails regions of the hydration
shell outside of interal cavities. A relatively robust workaround for
roughly globular proteins is to shrink the convex hull, using the
``scale`` argument of :func:`hop.qhull.make_ca_points`. Shrinking to
70% appears to be a good starting point::

  hop.qhull.points_from_selection(psf='protein.psf', pdb='protein.pdb', filename='ca_70%.dat', scale=0.7)

The convex hull itself is generated from the datafile of the points::

  Q70 = hop.qhull.ConvexHull('ca_70%.dat', workdir='cavity70%')

Another density grid ``b`` (such as a real water density for the bulk) is
currently required to generate a pseudo density based on the convex
hull. The real density provides the grid on which the convex hull is
mapped::

  b = hop.sitemap.Density(filename='bulk')
  QD70 = Q70.Density(b)

(This maps out sites at the threshold level set in ``b``; change it
with the :meth:`hop.sitemap.Density.map_sites` method if required.)

Insert a bulk density::
  QD70.site_insert_bulk(b)

"""

from __future__ import with_statement, absolute_import

import os
import errno
import tempfile
import shutil
from subprocess import Popen

import numpy

import logging

logger = logging.getLogger("MDAnalysis.app.qhull")

#: Comparisons of distances less than EPSILON yield equal.
EPSILON = 1e-6

def points_from_selection(*args, **kwargs):
    """Create a list of points from selected atoms in a format suitable for ``qhull``.

    points_from_selection(topology, structure, selection="name CA", filename="points.dat", scale=None)

    :Arguments:
    - psf: Charmm topology file
    - pdb: coordinates
    - selection: MDAnalysis select_atoms() selection string [C-alpha atoms]
    - filename: name of the output file; used as input for :class:`ConvexHull`
    - scale: scale points around the centre of geometry; values of 0.5 - 0.7 typically ensure that
      the convex hull is inside the protein; default is to not to scale, i.e. scale = 1.
    """

    from MDAnalysis import as_Universe
    u = as_Universe(*args, permissive=kwargs.pop('permissive', None))
    coordinates = u.select_atoms(kwargs.pop('selection', "name CA")).positions
    write_coordinates(kwargs.pop('filename', "points.dat"), coordinates, scale=kwargs.pop('scale',None))

def write_coordinates(filename, points, scale=None):
    """Write an array of points to a file suitable for qhull."""

    points = numpy.asarray(points)
    if not scale is None:
        center = points.mean(axis=0)
        points[:] = scale*points + (1-scale)*center  # scale*(points - center) + center
        logger.info("Scaled coordinates by factor %(scale)g relative to center of geometry %(center)r" % vars())

    with open(filename, 'w') as data:
        data.write('%d\n' % points.shape[1])  # dimension
        data.write('%d\n' % points.shape[0])  # number of points
        fmt = " ".join(["%f"]*points.shape[1]) + "\n"
        for point in points:
            data.write(fmt % tuple(point))
    logger.info("Wrote points to %(filename)r." % vars())


class ConvexHull(object):
    """The convex hull of a set of points.

    The convex hull is calculated with the `qhull`_ program.

    .. _qhull: http://www.qhull.org/
    """

    def __init__(self, coordinates, workdir=None, prefix=None):
        """Compute convex hull and populate data structures.

        :Arguments:
        - coordinates: input suitable for qconvex
        - workdir: store intermediate files in workdir (tmp dir by default)
        - prefix: filename prefix for intermediate output files
        """
        if workdir is None:
            self.workdir = tempfile.mkdtemp(prefix="tmp", suffix="_ConvexHull")
            self.workdir_is_temp = True
        else:
            self.workdir = workdir
            self.workdir_is_temp = False
            try:
                os.mkdir(self.workdir)
            except OSError, err:
                if err.errno != errno.EEXIST:
                    raise

        self.prefix = prefix or ""
        self.files = {'coordinates': coordinates,
                      'planes': self.wd(self.prefix+"planes.dat"),
                      'vertices': self.wd(self.prefix+"vertices.dat"),}

        # run qconvex
        args = ['n', 'TO', "'"+self.files['planes']+"'"]
        rc = self.qconvex(args)
        if rc != 0:
            raise OSError(rc, "qconvex failed computing planes, rc=%(rc)d" % vars(), self.files['planes'])
        self.planes = self.read_planes()
        logger.debug("Wrote %d planes to %r" % (len(self.planes), self.files['planes']))

        args = ['p', 'TO', "'"+self.files['vertices']+"'"]
        rc = self.qconvex(args)
        if rc != 0:
            raise OSError(rc, "qconvex failed computing vertices, rc=%(rc)d" % vars(), self.files['vertices'])
        self.vertices = self.read_vertices()
        logger.debug("Wrote %d vertices to %r" % (len(self.vertices), self.files['vertices']))

    def wd(self, *args):
        """Return path in workdir."""
        return os.path.join(self.workdir, *args)

    def qconvex(self, args):
        with open(self.files['coordinates']) as coord:
            # must use stdin, TI option cannot deal with dots in filenames, eg 'ca.dat'
            Q = Popen(['qconvex']+args, stdin=coord)
            rc = Q.wait()
        return rc

    def read_vertices(self):
        """Read vertices from qconvex p file.

        Numpy array of points [[x,y,z], ...]
        """
        return self._data_reader('vertices')

    def read_planes(self):
        """Read planes from qconvex n file.

        Numpy array [[n1,n2,n3,-p], ...] for planes n*x = -p.

        Planes are oriented and point outwards.
        """
        return self._data_reader('planes')

    def _data_reader(self, name):
        """Read simple data structures from qhull files.

        :Arguments:
        - name: keyword in self.files

        File format::
           dimension
           npoints
           x1 x2 x3 x4 ...
        """
        with open(self.files[name]) as data:
            dimension = int(data.readline())  # eg 3+1 for planes, 3 for vertices
            npoints = int(data.readline())
            a = []
            for line in data:
                a.append(map(float, line.strip().split()))
        if len(a) != npoints:
            raise IOError("Wrong number of datapoints %d, should be %d" % (len(a), npoints))
        return numpy.array(a)


    def point_inside(self, point):
        """Check if point [x,y,z] is inside the polyhedron defined by planes.

        Iff for all i: plane[i]([x,y,z]) = n*[x,y,z] + p < 0 <==> [x,y,z] inside

        (i.e. [x,y,z] is under *all* planes and the planes completely define the enclosed space
        """
        # crappy implementation, I am sure one can do this better with broadcasts
        # or a better algorithm, eg
        # http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
        # 1. shoot semi-inifinite ray. 2. count how many faces F it crosses
        # crosses iff (1) P under plane of F (2) projection of P on plane is inside F (2D problem)
        # (Problems if ray hits a vertex.)

        # check ALL faces.... :-p
        return numpy.all(numpy.dot(self.planes[:,:3], point) + self.planes[:,3] <= EPSILON)

    def points_inside(self, points):
        """Return bool array for all points:

        True: inside
        False: outside

        :Arguments:
        - points = [[x1,y1,z1], ...] or an iterator that supplies points
        - planes: normal forms of planes

        :Returns:
        Array with truth values such as [True, False, True, ...]
        """
        # crappy code, should optimize for numpy
        return numpy.fromiter((self.point_inside(point) for point in points), numpy.bool)

    def write_vertices_pdb(self, pdb="vertices.pdb"):
        ppw = VertexPDBWriter(pdb)
        ppw.write(self.vertices)
        logger.info("Wrote vertices to pdb file %(pdb)r." % vars())

    def Density(self, density, fillvalue=None):
        """Create a Density object of the interior of the convex hall.

        Uses another Density object *density*  as a template for the grid.

        .. Note:: This is rather slow and should be optimized.
        """

        from hop.sitemap import Density

        if fillvalue is None:
            fillvalue = 2*density.P['threshold']

        # TODO: OPTIMIZE
        # This is S-L-O-W because density.centers is slow (but at least a iterator using ndindex)
        # Reshaping relies on the centers being in correct order (provided by numpy.ndindex())
        mask = self.points_inside(density.centers()).reshape(density.grid.shape)

        grid = numpy.zeros_like(density.grid)
        grid[mask] = fillvalue        # fill the inside with high density
        parameters = density.P.copy()
        try:
            del parameters['bulk_site']
        except KeyError:
            pass
        return Density(grid=grid, edges=density.edges, parameters=parameters, unit=density.unit)

    def __del__(self):
        if self.workdir_is_temp:
            shutil.rmtree(self.workdir, ignore_errors=True)


class VertexPDBWriter(object):
    """PDB writer that implements a subset of the PDB 3.2 standard.
    http://www.wwpdb.org/documentation/format32/v3.2.html
    """
    #          1         2         3         4         5         6         7         8
    # 123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.
    # ATOM__seria nameAres CressI   xxxxxxxxyyyyyyyyzzzzzzzzOCCUPAtempft          elCH
    # ATOM  %5d   %-4s %-3s %4d %1s %8.3f   %8.3f   %8.3f   %6.2f %6.2f           %2s
    #                 %1s  %1s                                                      %2d
    #            =        =      ===                                    ==========
    # ATOM  %5d %-4s%1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2d
    # ATOM  %(serial)5d %(name)-4s%(altLoc)1s%(resName)-3s %(chainID)1s%(resSeq)4d%(iCode)1s   %(x)8.3f%(y)8.3f%(z)8.3f%(occupancy)6.2f%(tempFactor)6.2f          %(element)2s%(charge)2d

    fmt = {'ATOM':   "ATOM  %(serial)5d %(name)-4s%(altLoc)1s%(resName)-3s %(chainID)1s%(resSeq)4d%(iCode)1s   %(x)8.3f%(y)8.3f%(z)8.3f%(occupancy)6.2f%(tempFactor)6.2f          %(element)2s%(charge)2d\n",
           'REMARK': "REMARK     %s\n",
           'TITLE':  "TITLE    %s\n",
           }

    def __init__(self,filename):
        self.filename = filename
        self.pdb = open(self.filename,'w')

    def close(self):
        self.pdb.close()

    def write(self,coordinates,name="CA",resname="VRT",resid=1):
        """Write coordinates as CA."""

        self.TITLE("points as CA")
        for i, atom in enumerate(coordinates):
            self.ATOM(serial=i+1, name=name.strip(), resName=resname.strip(), resSeq=resid,
                      x=coordinates[i,0], y=coordinates[i,1], z=coordinates[i,2])
        self.close()

    def TITLE(self,*title):
        """Write TITLE record.
        http://www.wwpdb.org/documentation/format32/sect2.html
        """
        line = " ".join(title)    # should do continuation automatically
        self.pdb.write(self.fmt['TITLE'] % line)

    def REMARK(self,*remark):
        """Write generic REMARK record (without number).
        http://www.wwpdb.org/documentation/format32/remarks1.html
        http://www.wwpdb.org/documentation/format32/remarks2.html
        """
        line = " ".join(remark)
        self.pdb.write(self.fmt['REMARK'] % line)

    def ATOM(self,serial=None,name=None,altLoc=None,resName=None,chainID=None,
             resSeq=None,iCode=None,x=None,y=None,z=None,occupancy=1.0,tempFactor=0.0,
             element=None,charge=0):
        """Write ATOM record.
        http://www.wwpdb.org/documentation/format32/sect9.html
        Only some keword args are optional (altLoc, iCode, chainID), for some defaults are set.

        All inputs are cut to the maximum allowed length. For integer
        numbers the highest-value digits are chopped (so that the
        serial and reSeq wrap); for strings the trailing characters
        are chopped.

        Note: Floats are not checked and can potentially screw up the format.
        """
        for arg in ('serial','name','resName','resSeq','x','y','z',
                    'occupancy','tempFactor','charge'):
            if locals()[arg] is None:
                raise ValueError('parameter '+arg+' must be defined.')
        serial = int(str(serial)[-5:])  # check for overflow here?
        name = name[:4]
        if len(name) < 4:
            name = " "+name   # customary to start in column 14
        altLoc = altLoc or " "
        altLoc= altLoc[:1]
        resName = resName[:3]
        chainID = chainID or ""   # or should we provide a chainID such as 'A'?
        chainId = chainID[:1]
        resSeq = int(str(resSeq)[-4:]) # check for overflow here?
        iCode = iCode or ""
        iCode = iCode[:1]
        element = element or name.strip()[0]  # could have a proper dict here...
        element = element[:2]
        self.pdb.write(self.fmt['ATOM'] % vars())

    def __del__(self):
        self.close()

