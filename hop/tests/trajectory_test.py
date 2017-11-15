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
directly or file a enhancement request at http://github.com/orbeckst/hop/issues


Some common selection strings:

  * "name OW" for water in Gromacs
  * "name OH2" for water in CHARMM
"""

import os.path, errno
import numpy 
import MDAnalysis
import hop.interactive
import hop.trajectory
from hop.utilities import unlink_f, mkdir_p
from six.moves import zip

import logging
logger = logging.getLogger('MDAnalysis.app')


if __name__ == "__main__":
    import sys
    from optparse import OptionParser

    parser = OptionParser(usage=__doc__)
    parser.add_option("-f", "--trajdcd", dest="trajdcd",
                      metavar="DCD",
                      help="Trajectory dcd file [%default]")
    parser.add_option('-s', '--trajpsf', dest="trajpsf",
                     metavar="PSF",
                     help="Trajectory psf file [%default]")
    parser.add_option("-D", "--analysisdir", dest="analysisdir",
                      metavar="DIRNAME",
                      help="results will be stored under DIRNAME/(basename DIR)  [%default]")
    parser.add_option("-t", "--analysistrajdcd", dest="analysistrajdcd",
                      metavar="ADCD",
                      help="Comparison trajectory dcd file  [%default]")
    parser.add_option("-y", "--analysistrajpsf", dest="analysistrajpsf",
                      metavar="APSF",
                      help="Comparison trajectory psf file  [%default]")

    opts,args = parser.parse_args()

    analysisdir = os.path.join(opts.analysisdir)
    
    #new_trajectory = hop.trajectory.HoppingTrajectory(hopdcd=opts.trajdcd,hoppsf=opts.trajpsf)
    #old_trajectory = hop.trajectory.HoppingTrajectory(filename=opts.analysistrajectory)
    
    new_trajectory = MDAnalysis.Universe(opts.trajpsf,opts.trajdcd)
    old_trajectory = MDAnalysis.Universe(opts.analysistrajpsf,opts.analysistrajdcd)

    MDAnalysis.start_logging()

    if len(args) != 0:
        logger.fatal("This command only accepts option arguments. See --help.")
        sys.exit(1)

    trajectorydcd = os.path.abspath(opts.trajdcd)
    if not os.path.exists(trajectorydcd):
        errmsg = "Trajectory dcd %(trajectorydcd)r not found; (use --trajdcd)" % vars()
        logger.fatal(errmsg)
        raise IOError(errno.ENOENT, errmsg)

    trajectorypsf = os.path.abspath(opts.trajpsf)
    if not os.path.exists(trajectorypsf):
        errmsg = "Trajectory psf %(trajectorypsf)r not found; (use --trajpsf)" % vars()
        logger.fatal(errmsg)
        raise IOError(errno.ENOENT, errmsg)

    analysistrajdcd = os.path.abspath(opts.analysistrajdcd)
    if not os.path.exists(analysistrajdcd):
        errmsg = "Analysis trajectory dcd %(analysistrajdcd)r not found; (use --analysistraj)" % vars()
        logger.fatal(errmsg)
        raise IOError(errno.ENOENT, errmsg)

    analysistrajpsf = os.path.abspath(opts.analysistrajpsf)
    if not os.path.exists(analysistrajpsf):
        errmsg = "Analysis trajectory psf %(analysistrajpsf)r not found; (use --analysistraj)" % vars()
        logger.fatal(errmsg)
        raise IOError(errno.ENOENT, errmsg)
    #logger.debug("trajectory   = %(trajectory)r", vars())

    startdirectory = os.path.abspath(os.path.curdir)
     
    ref = old_trajectory
    u = new_trajectory
    
    for ts_reference, ts_test in zip(ref.trajectory, u.trajectory):
        numpy.testing.assert_array_almost_equal(ts_test.positions[:, 0], ts_test.positions[:, 0],
                   err_msg="Hop trajectory sites do not agree at frame {0}".format(ts_reference.frame)) 
MDAnalysis.stop_logging()

