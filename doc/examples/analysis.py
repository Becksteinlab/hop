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
job = Job(inputfiles=dict(dens_apo  = '1IFC/analysis/water.pickle',
                          bulk_apo =  '1IFC/analysis/bulk.pickle',
                          dens_holo = '2IFB/analysis/water.pickle',
                          bulk_holo = '2IFB/analysis/bulk.pickle',
                          ),
          outputfiles=dict(scan = 'analysis/scan.pickle',
                           scan_pdf = 'figs/scan.pdf',
                           ))
#
#------------------------------------------------------------
job.stage()
F = job.filenames

# commands
import hop.utilities    # must come first (for reasons unknown)
hop.utilities.matplotlib_interactive(False)  # no X11 available

import hop.analysis
import numpy

DA =  hop.analysis.DensityAnalysis(F['dens_holo'],reference=F['dens_apo'],
                                   bulkdensities=F['bulk_holo'],
                                   refbulkdensity=F['bulk_apo'])
DA.scan(numpy.arange(0.6,4.5,0.05))
DA.save(F['scan'])
DA.plotscan(F['scan_pdf'])

job.unstage()
job.cleanup()
