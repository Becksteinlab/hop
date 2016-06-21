#!/usr/bin/env python
#---------------- EDIT JOB NAME -------------------
#$ -N hoptraj
#--------------------------------------------------
#$ -S /usr/bin/python
#$ -v PYTHONPATH=/home/oliver/Library/python-lib
#$ -v LD_LIBRARY_PATH=/opt/intel/cmkl/8.0/lib/32:/opt/intel/itc60/slib:/opt/intel/ipp41/ia32_itanium/sharedlib:/opt/intel/ipp41/ia32_itanium/sharedlib/linux32:/opt/intel/fc/9.0/lib:/opt/intel/cc/9.0/lib
#$ -r n
#$ -j y
# Using the current working directory is IMPORTANT with the default settings for Job()
#$ -cwd
#$ -m e
# $Id$

        
from staging.SunGridEngine import Job

#------------------------------------------------------------
# EDIT THE inputfiles AND outputfiles DICTIONARIES.
#------------------------------------------------------------
# record input and output files relative to top_dir = cwd
job = Job(inputfiles=dict(density='analysis/water.pickle',
                          bulk='analysis/bulk.pickle',
                          psf='inp/XXX.psf',
                          dcd='trj/rmsfit_XXX.dcd'),
          outputfiles=dict(density='analysis/water.pickle',
                           hopdcd='trj/hoptraj.dcd',
                           hoppsf='trj/hoptraj.psf',
                           ))
#
#------------------------------------------------------------

job.stage()

# commands
import hop.utilities
hop.utilities.matplotlib_interactive(False)
from hop.interactive import *
from hop.sitemap import Density
import hop.constants

density = Density(filename=job.input['density'])
bulk = Density(filename=job.input['bulk'])

# Fixing metadata -- only necessary if the dcd has a wrong header
# as for instance produced by catdcd.
##density.metadata['psf'] = job.input['psf'] # used in make_hoppingtraj()
##density.metadata['dcd'] = job.input['dcd']
##delta_ps = 0.5   # the time between two frames in ps (Hirsh's trajectories)
##delta_AKMA = delta_ps * hop.constants.get_conversion_factor('time','ps','AKMA')
##density.metadata['delta'] = delta_AKMA
##fixtrajectory = {'delta':density.metadata['delta']}
fixtrajectory = None

# Add the biggest bulk site at position 1 if we haven't done so already.
# This is important so we are making extra sure.
try:
    density.site_insert_bulk(bulk)   # hack!
except ValueError,errmsg:
    print errmsg
    print "Bulk site not inserted because there already exists a bulk site --- that's good!"
density.save(job.output['density'])  # also save modified metadata
del bulk

hops = make_hoppingtraj(density,job.output['hopdcd'],fixtrajectory=fixtrajectory)

job.unstage()
job.cleanup()
