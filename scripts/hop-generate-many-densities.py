#!/usr/bin/env python
"""%prog [options] DIR [DIR ...]

Generate densities (solvent and bulk) for identically laid out simulations in
DIR. The assumption is that each simulation lives in a separate directory and
that the directory name distinguishes simulations. All simulations have the
same topology and the trajectory name is the same in each DIR.

At the moment, default values are used for most settings (because this is a
primitive script). In particular:

 * The output files are always named "water.pickle" and "bulk.pickle" and they
   are stored under the analysis dir.
 * All options for setting up the bulk are fixed.
 * The grid spacing is fixed at 1 A.

For more fine grained control, use hop.interactive.generate_densities()
directly.

Some common selection strings:

  * "name OW" for water in Gromacs
  * "name OH2" for water in CHARMM
"""

import os.path, errno
import MDAnalysis
import hop.interactive
from hop.utilities import unlink_f, mkdir_p

import logging
logger = logging.getLogger('MDAnalysis.app.hop')



def generate_densities_locally(topology, trajectory, atomselection, localcopy=False):
    def _generate_densities(traj):
        return hop.interactive.generate_densities(topology, traj, atomselection=atomselection)
    if localcopy:
        from tempfile import mkstemp
        from shutil import copy
        root,ext = os.path.splitext(trajectory)
        fd, tmptrajectory = mkstemp(suffix=ext)
        logger.info("Making local copy to improve read performance: %(trajectory)r --> %(tmptrajectory)r" % vars())
        try:
            copy(trajectory, tmptrajectory)
            densities = _generate_densities(tmptrajectory)
        finally:
            unlink_f(tmptrajectory)
    else:
        densities = _generate_densities(trajectory)
    return densities

if __name__ == "__main__":
    import sys
    from optparse import OptionParser

    parser = OptionParser(usage=__doc__)
    parser.add_option("-s", "--topology", dest="topology",
                      metavar="FILE",
                      help="topology to go with the trajectories; can be a PSF "
                      "PDB, GRO, or whatever else MDAnalysis accepts [%default]")
    parser.add_option("-f", "--trajectory", dest="trajectory",
                      metavar="FILE",
                      help="rms-fitted trajectory name to be found under DIR [%default]")
    parser.add_option("-A", "--atom-selection", dest="atomselection",
                      metavar="SELECTION",
                      help="MDAnalysis selection string to pick which atoms are being counted; "
                      "for water one typically chooses the water oxygen. The value depends crucially "
                      "on the atom names defined in the topology [%default].")
    parser.add_option("-D", "--analysisdir", dest="analysisdir",
                      metavar="DIRNAME",
                      help="results will be stored under DIRNAME/(basename DIR)  [%default]")
    parser.add_option("-l", "--local-copy", dest="localcopy",
                      action='store_true',
                      help="copy trajectory to a temporary local disk for better read performance. "
                      "Requires sufficient space in TEMP.")

    parser.set_defaults(topology="md.pdb", trajectory="rmsfit_md.xtc",
                        atomselection="name OW",
                        analysisdir="analysis")

    opts,args = parser.parse_args()

    MDAnalysis.start_logging()

    if len(args) == 0:
        logger.fatal("At least one directory is required. See --help.")
        sys.exit(1)

    topology = os.path.abspath(opts.topology)
    if not os.path.exists(topology):
        errmsg = "Topology %(topology)r not found; (use --topology)" % vars()
        logger.fatal(errmsg)
        raise IOError(errno.ENOENT, errmsg)

    startdirectory = os.path.abspath(os.curdir)
    for d in args:
        os.chdir(startdirectory)
        logger.info("Generating densities for dir %(d)r", vars())
        if not os.path.exists(d):
            logger.fatal("Directory %r does not exist.", d)
            raise IOError(errno.ENOENT, d)
        analysisdir = os.path.join(opts.analysisdir, os.path.basename(d))
        try:
            mkdir_p(analysisdir)
        except:
            logger.exception()
            raise
        trajectory = os.path.abspath(os.path.join(d, opts.trajectory))
        if not os.path.exists(trajectory):
            errmsg = "Trajectory %(trajectory)r not found; (use --trajectory)" % vars()
            logger.fatal(errmsg)
            raise IOError(errno.ENOENT, errmsg)
        logger.debug("topology   = %(topology)r", vars())
        logger.debug("trajectory = %(trajectory)r", vars())
        logger.debug("selection  = %(atomselection)r", vars(opts))

        try:
            os.chdir(analysisdir)
            logger.debug("Working in %(analysisdir)r..." % vars())
            densities = generate_densities_locally(topology, trajectory, opts.atomselection,
                                                   localcopy=opts.localcopy)
        finally:
            os.chdir(startdirectory)

    MDAnalysis.stop_logging()

