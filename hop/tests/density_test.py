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
import hop.density
from hop.utilities import unlink_f, mkdir_p

import logging
logger = logging.getLogger('MDAnalysis.app')


if __name__ == "__main__":
    import sys
    from optparse import OptionParser

    parser = OptionParser(usage=__doc__)
    parser.add_option("-s", "--density", dest="density",
                      metavar="FILE",
                      help="Density file (e.g. water.pickle) for analysis [%default]")
    parser.add_option("-D", "--analysisdir", dest="analysisdir",
                      metavar="DIRNAME",
                      help="results will be stored under DIRNAME/(basename DIR)  [%default]")
    parser.add_option("-t", "--analysisdensity", dest="analysisdensity",
                      metavar="DIRNAME",
                      help="Comparison density file  [%default]")

    opts,args = parser.parse_args()

    analysisdir = os.path.join(opts.analysisdir)
    
    new_density = hop.sitemap.Density(filename=opts.density)
    old_density = hop.sitemap.Density(filename=opts.analysisdensity)

    MDAnalysis.start_logging()

    if len(args) != 0:
        logger.fatal("This command only accepts option arguments. See --help.")
        sys.exit(1)

    density = os.path.abspath(opts.density)
    if not os.path.exists(density):
        errmsg = "Density %(density)r not found; (use --density)" % vars()
        logger.fatal(errmsg)
        raise IOError(errno.ENOENT, errmsg)

    analysisdensity = os.path.abspath(opts.analysisdensity)
    if not os.path.exists(analysisdensity):
        errmsg = "Analysis density %(analysisdensity)r not found; (use -t)" % vars()
        logger.fatal(errmsg)
        raise IOError(errno.ENOENT, errmsg)

    logger.debug("density   = %(density)r", vars())

    startdirectory = os.path.abspath(os.path.curdir)
     
    center_original=[i for i in old_density.centers()]

    center_new=[i for i in new_density.centers()]

    numpy.testing.assert_array_almost_equal(center_original,center_new,err_msg="Centers are unequal")

    del center_original, center_new

    
    
#    numpy.testing.assert_array_almost_equal(new_density.site_properties,old_density.site_properties,err_msg="Site Properties unequal")
    '''
    if new_density.site_properties != old_density.site_properties:
        print('Site properties are unequal')

    ''' 
    if new_density.sites != old_density.sites:
        print("Site lists Unequal")
    #if numpy.testing.assert_array_almost_equal(new_density.sites,old_density.sites) is False:
    #    print('different sites')


    numpy.testing.assert_array_almost_equal(new_density.site_volume()[1],numpy.sort(new_density.site_volume()[1])[::-1],err_msg='Sites not sorted') 

    numpy.testing.assert_array_almost_equal(new_density.site_volume(),old_density.site_volume(),err_msg='Site volumes are unequal') 

    numpy.testing.assert_array_almost_equal(new_density.site_occupancy(),old_density.site_occupancy(),err_msg="Site occupancies are unequal")

    numpy.testing.assert_array_almost_equal(new_density.map,old_density.map,err_msg="Maps are unequal")

    if numpy.testing.assert_array_almost_equal(new_density.origin,old_density.origin) is False:
        print('different origins')
MDAnalysis.stop_logging()

