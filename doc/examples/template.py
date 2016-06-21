#!/usr/bin/env python
#---------------- EDIT JOB NAME -------------------
#$ -N WaterAnalysis
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
job = Job(variables=dict(),
          inputfiles=dict(
                          ),
          outputfiles=dict(
                           ))
#
#------------------------------------------------------------
job.stage()
F = job.filenames
V = job.variables

# commands
import hop.utilities    # must come first (for reasons unknown)
hop.utilities.matplotlib_interactive(False)  # no X11 available

import hop.analysis



job.unstage()
job.cleanup()
