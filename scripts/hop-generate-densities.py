#!/usr/bin/env python
"""%prog [options] -s TOPOL -f TRAJECTORY

Generate densities (solvent and bulk) for the system specified by
the structure file TOPOL and the MD TRAJECTORY. Any combination of
TOPOL and TRAJECTORY that can be read by MDAnalysis is acceptable
(e.g. a PSF/DCD or GRO/XTC combination).

At the moment, default values are used for most settings (because this is a
primitive script). In particular:

 * The output files are always named "water.pickle" and "bulk.pickle" and they
   are stored under the analysis dir.
 * The density threshold for defining bulk is fixed at exp(-0.5) = 0.60...

For more fine grained control, use hop.interactive.generate_densities()
directly or file a enhancement request at http://github.com/Becksteinlab/hop/issues


Some common selection strings:

  * "name OW" for water in Gromacs
  * "name OH2" for water in CHARMM
"""

import os.path, errno
import numpy
import MDAnalysis
import hop.interactive
from hop.utilities import unlink_f, mkdir_p

import logging
logger = logging.getLogger('MDAnalysis.app')



def generate_densities_locally(topology, trajectory, localcopy=False, **kwargs):
    # default values in hop.density.DensityCreator
    def _generate_densities(traj):
        return hop.interactive.generate_densities(topology, traj, **kwargs)
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
    parser.add_option("--threshold", type='float',dest="solvent_threshold",
                      metavar="CUTOFF",
                      help="hydration sites are considered regions of density greater than CUTOFF; "
                      "CUTOFF is measured in units of the water density at standard confitions; "
                      "e.g. 2.0 means twice the bulk water density. [%default]")
    parser.add_option("--delta", dest="delta",
                      metavar="DELTA",
                      help="The density is histogrammed on a grid of spacing DELTA Angstroem [%default]")
    parser.add_option("--bulk-solvent-distance", dest="cutoff",
                      metavar="CUTOFF",
                      help="bulk-water is assumed to start at a distance CUTOFF Angstroem from the "
                      "solute selection SELECTION [%default]")
    parser.add_option("--bulk-solute-selection", dest="soluteselection",
                      metavar="SELECTION",
                      help="MDAnalysis selection string to select the solute (for bulk "
                      "density)  [%default]")
    parser.add_option("--local-copy", dest="localcopy",
                      action='store_true',
                      help="copy trajectory to a temporary local disk for better read performance. "
                      "Requires sufficient space in TEMP.")

    parser.set_defaults(topology="md.pdb", trajectory="rmsfit_md.xtc",
                        atomselection="name OW", solvent_threshold=numpy.e, delta=1.0,
                        cutoff=3.5, soluteselection="protein and not name H*",
                        analysisdir="analysis")

    opts,args = parser.parse_args()

    MDAnalysis.start_logging()

    if len(args) != 0:
        logger.fatal("This command only accepts option arguments. See --help.")
        sys.exit(1)

    topology = os.path.abspath(opts.topology)
    if not os.path.exists(topology):
        errmsg = "Topology %(topology)r not found; (use --topology)" % vars()
        logger.fatal(errmsg)
        raise IOError(errno.ENOENT, errmsg)
    trajectory = os.path.abspath(opts.trajectory)
    if not os.path.exists(trajectory):
        errmsg = "Trajectory %(trajectory)r not found; (use --trajectory)" % vars()
        logger.fatal(errmsg)
        raise IOError(errno.ENOENT, errmsg)

    analysisdir = os.path.join(opts.analysisdir)
    try:
        mkdir_p(analysisdir)
    except:
        logger.exception()
        raise
    logger.debug("topology   = %(topology)r", vars())
    logger.debug("trajectory = %(trajectory)r", vars())

    for v in ('atomselection', 'solvent_threshold', 'delta', 'cutoff', 'soluteselection'):
        fmt = "%(v)-17s = %%(%(v)s)s" % vars()
        logger.debug(fmt, vars(opts))

    startdirectory = os.path.abspath(os.path.curdir)
    try:
        os.chdir(analysisdir)
        logger.debug("Working in %(analysisdir)r..." % vars())
        densities = generate_densities_locally(topology, trajectory, atomselection=opts.atomselection,
                                               solvent_threshold=opts.solvent_threshold, delta=opts.delta,
                                               soluteselection=opts.soluteselection, cutoff=opts.cutoff,
                                               localcopy=opts.localcopy)
    finally:
        os.chdir(startdirectory)

    MDAnalysis.stop_logging()

