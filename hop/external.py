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

__doc__ = """Generate densities with the help of external applications.

Currently, only VMD is implemented via vmd.control. We launch VMD in
text mode and start a simple server, then send commands to the server
to run the volmap plugin, and finally pick up the dx file.
"""

import os,os.path, errno
import vmd

class VMD:
    """Launch VMD in text mode and start a simple server, then send
    commands to the server to run the volmap plugin, and finally pick up
    the dx file."""
    
    def __init__(self,**kwargs):
        """Start the VMD server process."""
        self.server = vmd.control.server(**kwargs)  # start server
        self.command = self.server.command

    def volmap(self, psf,dcd,dx, delta=1.0,load_new=True,
               atomselection='name OH2', **sel_args):
        """Calculate the water density around the protein with VMD's VolMap.

        volmap(psf,dcd,dx,delta=1.0,atomselection=<VMD atomselect>,**kwargs)        

        Arguments:
        
        psf       psf topology file
        dcd       trajectory that is RMS-fitted on the protein
        dx        output file for the density (OpenDX format), in A^-3.
        delta     size of a grid cell in Angstrom
        load_new  True: load dx and psf in VMD, False: just send volmap command
        atomselection
                  VMD selection string (with python keyword-interpolation of
                  all additional sel_args)

        This function just loads the data into VMD and runs the VolMap
        plugin with appropriate values
        (http://www.ks.uiuc.edu/Research/vmd/current/ug/node141.html and
        http://www.ks.uiuc.edu/Research/vmd/plugins/volmapgui/).

        The default atom selection being used calculates the water density in A^-3:
           atomselect top {name OH2}

        The VolMap checkpoint:*.dx file is automatically removed.
        
        Examples:

        * Bulk water density (exclude water near protein)

            atomselection='name OH2 and not within %(exclusion)f of protein',
            exclusion=3.5
        """

        # build selection string
        atomselectionstring = atomselection % sel_args
        
        # must convert paths to absolute paths as we don't know where the
        # server's cwd is
        psf,dcd,dx = map(os.path.abspath, [psf,dcd,dx])

        if load_new:
            c = self.command(\
               'mol new %(psf)s type psf first 0 last -1 step 1 filebonds 1 '
               'autobonds 0 waitfor all' % locals(),
               'mol addfile %(dcd)s type dcd first 0 last -1 step 1 filebonds 1 '
               'autobonds 0 waitfor all' % locals())
        c = self.command(\
            'volmap density [atomselect top {%(atomselectionstring)s}] '
            '-res %(delta)f -allframes -combine avg -radscale 1.0 '
            '-o %(dx)s' % locals())

        # clean up
        path,filename = os.path.split(dx)
        checkpoint = os.path.join(path,'checkpoint:'+filename)
        try:
            os.remove(checkpoint)
        except os.error,e:
            if e.errno != errno.ENOENT:
                raise
        
        return c.results()
        

