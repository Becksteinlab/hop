#!/usr/bin/env python
#---------------- EDIT JOB NAME -------------------
#$ -N MCMC
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
job = Job(variables=dict(state='apo'),
          inputfiles=dict(hopgraph='analysis/hopgraph.pickle'
                          ),
          outputfiles=dict(scan = 'analysis/pscan.pickle',
                           occupancy_pdf = 'figs/mcmc_occupancy.pdf',
                           correl_pdf = 'figs/mcmc_occupancy_correl.pdf',                           
                           ))
#
#------------------------------------------------------------
job.stage()
F = job.filenames
V = job.variables

# commands
import hop.utilities    # must come first (for reasons unknown)
hop.utilities.matplotlib_interactive(False)  # no X11 available

import hop.MCMC
from pylab import *

M = hop.MCMC.Pscan(Ntotal=1e7)
M.save(F['scan'])

figure(1)
M.plot_occupancy()
title('occupancy (I-FABP %(state)s)' % V)
savefig(F['occupancy_pdf'])

figure(2)
M.plot_correl()
title('occupancy correlation with MD (I-FABP %(state)s)' % V)
xlim(0,1.05)
ylim(0,1.02)
savefig(F['correl_pdf'])


job.unstage()
job.cleanup()
