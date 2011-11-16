#!/usr/bin/env python
"""%prog [options] -s HOPTRAJ.psf -f HOPTRAJ.dcd DENSITY

Generate the hopping graph from HOPTRAJ.psf and HOPTRAJ.dcd and
DENSITY. DENSITY must be the pickle file used to make the hopping
trajectory.

"""

import os.path, errno
import MDAnalysis
import hop.sitemap, hop.trajectory, hop.interactive
from hop.utilities import unlink_f, mkdir_p

import logging
logger = logging.getLogger('MDAnalysis.app')

def generate_hopgraph(topology, trajectory, densityfile, filename, localcopy=False, **hopargs):
    density = hop.sitemap.Density(filename=densityfile)
    def _generate_hopgraph(trajectory):
        hoptraj = hop.trajectory.HoppingTrajectory(hoppsf=topology, hopdcd=trajectory)
        tgraph = hop.interactive.build_hoppinggraph(hoptraj,density)
        return tgraph
    if localcopy:
        # should probably also write locally and then copy back
        from tempfile import mkstemp
        from shutil import copy
        root,ext = os.path.splitext(trajectory)
        fd, tmptrajectory = mkstemp(suffix=ext)
        logger.info("Making local copy to improve read performance: %(trajectory)r --> %(tmptrajectory)r" % vars())
        try:
            copy(trajectory, tmptrajectory)
            tgraph = _generate_hopgraph(tmptrajectory)
        finally:
            unlink_f(tmptrajectory)
    else:
        tgraph = _generate_hopgraph(trajectory)

    h = tgraph.hopgraph       # main result is the 'hopgraph'
    h.save(filename)          # save result
    logger.info("Saved hopgraph as %(filename)s.pickle", vars())

    hop.interactive.hopgraph_basic_analysis(h, density, filename)

    return h

def analyze_hopgraph(densityfile, filename):
    density = hop.sitemap.Density(filename=densityfile)
    h = hop.graph.HoppingGraph(filename=filename)
    hop.interactive.hopgraph_basic_analysis(h, density, filename)
    return h


if __name__ == "__main__":
    import sys
    from optparse import OptionParser

    parser = OptionParser(usage=__doc__)
    parser.add_option("-s", "--topology", dest="topology",
                      metavar="FILE",
                      help="PSF topology to go with the hopping trajectories")
    parser.add_option("-f", "--trajectory", dest="trajectory",
                      metavar="FILE",
                      help="hopping trajectory DCD")
    parser.add_option("-o", "--output-name", dest="output",
                      metavar="FILENAME",
                      help="the hopping graph FILENAME.pickle etc...; its path component is assumed "
                      "to be an analysis directory where we can write other files [%default]")
    parser.add_option("-l", "--local-copy", dest="localcopy",
                      action='store_true',
                      help="copy trajectory to a temporary local disk for better read performance. "
                      "Requires sufficient space in TEMP.")
    parser.add_option("--force", dest="force", action="store_true",
                      help="force rerunning graph generation; by default only the basic analysis is run "
                      "if the hopgraph file aleady exists")

    parser.set_defaults(topology='hoptraj.psf', trajectory='hoptraj.dcd', output="hopgraph")

    opts,args = parser.parse_args()

    MDAnalysis.start_logging()

    if len(args) == 0:
        logger.fatal("A pickled density with bulk site is required. See --help.")
        sys.exit(1)

    density = os.path.abspath(args[0])
    if not os.path.exists(density):
        errmsg = "Density %(density)r not found" % vars()
        logger.fatal(errmsg)
        raise IOError(errno.ENOENT, errmsg)

    # I think this is how we construct filenames in hop... see utilities.filename_function
    hopgraph_filename = os.path.splitext(opts.output)[0] + ".pickle"
    if not os.path.exists(hopgraph_filename) or opts.force:
        topology = os.path.abspath(opts.topology)
        if not os.path.exists(topology):
            errmsg = "Hopping topology %(topology)r not found; (use --topology)" % vars()
            logger.fatal(errmsg)
            raise IOError(errno.ENOENT, errmsg)
        trajectory = os.path.abspath(opts.trajectory)
        if not os.path.exists(trajectory):
            errmsg = "Hopping trajectory %(trajectory)r not found; (use --trajectory)" % vars()
            logger.fatal(errmsg)
            raise IOError(errno.ENOENT, errmsg)

        logger.info("Generating hopping graph from trajectory %(topology)r/%(trajectory)r ...", vars())
        logger.debug("density    = %(density)r", vars())
        logger.debug("topology   = %(topology)r", vars())
        logger.debug("trajectory = %(trajectory)r", vars())

        hopgraph = generate_hopgraph(topology, trajectory, density, opts.output,
                                     localcopy=opts.localcopy)
        logger.info("Created hopping graph %(output)s.pickle and other files", vars(opts))
    else:
        logger.info("hopgraph file %(hopgraph_filename)r already exists; only doing analysis", vars())
        analyze_hopgraph(density, opts.output)

    MDAnalysis.stop_logging()

