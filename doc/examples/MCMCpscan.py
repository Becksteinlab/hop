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
# $Id: markovsampling.py 2044 2008-07-22 19:21:12Z oliver $
"""Run 20 independent Pscans to check convergence of MCMC model with Ntotal."""
        
from staging.SunGridEngine import Job

#------------------------------------------------------------
# EDIT THE inputfiles AND outputfiles DICTIONARIES.
#------------------------------------------------------------
# record input and output files relative to top_dir = cwd
job = Job(variables=dict(repeats=20,
                         parameter='Ntotal',
                         pvalues=[    3e3,5e3,8e3, 
                                  1e4,3e4,5e4,8e4,
                                  1e5,3e5,5e5,8e5,
                                  1e6,5e6,
                                  1e7,3e7],
                         ),
          inputfiles=dict(hopgraph='analysis/hopgraph.pickle'
                          ),
          outputfiles=dict(scan = 'analysis/pscan_Ntotal.pickle',
                           scanfig = 'figs/pscan_Ntotal.pdf',
                           )
          )
#
#------------------------------------------------------------
job.stage()
F = job.filenames
V = job.variables

# commands
import hop.utilities    # must come first (for reasons unknown)
hop.utilities.matplotlib_interactive(False)  # no X11 available

import hop.MCMC
import pylab

mpscan = hop.MCMC.MultiPscan(repeat=V['repeats'],
                             parameter=V['parameter'],pvalues=V['pvalues'],
                             filename=F['hopgraph'])
mpscan.save(F['scan'])
mpscan.plot()

pylab.savefig(F['scanfig'])

job.unstage()
job.cleanup()

