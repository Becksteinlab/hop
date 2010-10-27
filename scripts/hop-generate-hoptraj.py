#!/usr/bin/env python
"""%prog [options] DENSITY

Generate the hopping trajectory from the DENSITY and the original
trajectory. DENSITY must be a pickle file produced by hop and it
should contain the bulk site; otherwise most other downstream analysis
makes little sense.

If ``-s topology`` and ``-f trajectory`` are left out then the values
for *topology* and *trajectory* stored in DENSITY are used; the
paths might be incorrect relative to the starting directory so it is
recommended to supply them on the command line.
"""

import os.path, errno
import MDAnalysis
import hop.sitemap, hop.trajectory
from hop.utilities import unlink_f, mkdir_p

import logging
logger = logging.getLogger('MDAnalysis.app')



def generate_hoptraj_locally(topology, trajectory, density, filename, atomselection,
                             localcopy=False, **hopargs):
    def _generate_hoptraj(traj):
        try:
            if len(density.sites) < 2:
                raise ValueError
        except AttributeError,ValueError:
            errmsg = 'The density misses a site map or has only one site.'
            logger.fatal(errmsg)
            raise ValueError(errmsg)
        u = MDAnalysis.Universe(topology, traj)
        group = u.selectAtoms(atomselection)
        hops = hop.trajectory.HoppingTrajectory(u.trajectory,group,density,**hopargs)
        hops.write(filename)
        return hops

    if localcopy:
        from tempfile import mkstemp
        from shutil import copy
        root,ext = os.path.splitext(trajectory)
        fd, tmptrajectory = mkstemp(suffix=ext)
        logger.info("Making local copy to improve read performance: %(trajectory)r --> %(tmptrajectory)r" % vars())
        try:
            copy(trajectory, tmptrajectory)
            hops = _generate_hoptraj(tmptrajectory)
        finally:
            unlink_f(tmptrajectory)
    else:
        hops = _generate_hoptraj(trajectory)
    return hops

if __name__ == "__main__":
    import sys
    from optparse import OptionParser

    parser = OptionParser(usage=__doc__)
    parser.add_option("-s", "--topology", dest="topology",
                      metavar="FILE",
                      help="topology to go with the trajectories; can be a PSF "
                      "PDB, GRO, or whatever else MDAnalysis accepts; default is to "
                      "try the path recorded in DENSITY")
    parser.add_option("-f", "--trajectory", dest="trajectory",
                      metavar="FILE",
                      help="rms-fitted trajectory (default is to take the path recorded in DENSITY)")
    parser.add_option("-A", "--atom-selection", dest="atomselection",
                      metavar="SELECTION",
                      help="MDAnalysis selection string to pick which atoms are being counted; "
                      "for water one typically chooses the water oxygen. The value depends crucially "
                      "on the atom names defined in the topology. The default is to use the value "
                      "recorded in DENSITY.")
    parser.add_option("-o", "--output-name", dest="output",
                      metavar="FILENAME",
                      help="the hopping trajectory FILENAME.dcd and FILENAME.psf [%default]")
    parser.add_option("-l", "--local-copy", dest="localcopy",
                      action='store_true',
                      help="copy trajectory to a temporary local disk for better read performance. "
                      "Requires sufficient space in TEMP.")

    parser.set_defaults(topology=None, trajectory=None, 
                        atomselection=None,
                        output="hoptraj")

    opts,args = parser.parse_args()

    if len(args) == 0:
        logger.fatal("A pickled density with bulk site is required. See --help.")
        sys.exit(1)

    density = hop.sitemap.Density(filename=args[0])

    MDAnalysis.start_logging()
    if opts.topology:
        topology = os.path.abspath(opts.topology)
    else:
        topology = os.path.abspath(density.metadata['PSF'])
    if not os.path.exists(topology):
        errmsg = "Topology %(topology)r not found; (use --topology)" % vars()
        logger.fatal(errmsg)
        raise IOError(errno.ENOENT, errmsg)

    if opts.trajectory:
        trajectory = os.path.abspath(opts.trajectory)
    else:
        trajectory = os.path.abspath(density.metadata['DCD'])
    if not os.path.exists(trajectory):
        errmsg = "Trajectory %(trajectory)r not found; (use --trajectory)" % vars()
        logger.fatal(errmsg)
        raise IOError(errno.ENOENT, errmsg)

    if opts.atomselection:
        atomselection = opts.atomselection
    else:
        atomselection = density.metadata['atomselection']

    #startdirectory = os.path.abspath(os.curdir)
    #os.chdir(startdirectory)            

    logger.info("Generating hopping trajectory for density %r", args[0])
    logger.debug("density    = %r", args[0])
    logger.debug("topology   = %(topology)r", vars())
    logger.debug("trajectory = %(trajectory)r", vars())
    logger.debug("selection  = %(atomselection)r", vars()) 

    if not density.has_bulk():
        raise ValueError("The density does not have a bulk site---insert one!")

    hoptraj = generate_hoptraj_locally(topology, trajectory, density, opts.output,
                                         atomselection, localcopy=opts.localcopy)
    logger.info("Created hopping trajectory %(output)s.dcd with %(output)s.psf", vars(opts))
    MDAnalysis.stop_logging()        

