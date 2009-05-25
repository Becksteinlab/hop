#!/usr/bin/env python
#---------------- EDIT JOB NAME -------------------
#$ -N rmsfit1OPA
#--------------------------------------------------
#$ -S /usr/bin/python
# Importing the whole environment is important for setting up Charmm:
#$ -V
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
job = Job(inputfiles=dict(HEADER='str/header.str',
                          PSF='inp/crbp_apo.psf',
                          DCD='trj/1opa_salt_ewald_shake_10ang_prod.dcd',    
                          REF_PDB='coord/rmsfit_1opa_a_xtal.pdb',
                          ),
          outputfiles=dict(DCD_RMS='trj/rmsfit_1opa_salt_ewald_shake_10ang_prod.dcd',
                           CHARMM_LOG='log/rmsfit.log',
                           ))
#
#------------------------------------------------------------

job.stage()

# commands
import os

# This setup requires my canonical filesystem layout under CHARMM:
if not 'CHARMM' in os.environ:   # get sane value for CHARMM top dir
    os.environ['CHARMM'] = os.path.join(os.environ['BIO_D'],'Library','Charmm')
    os.environ['CHARMM_TOP'] = os.path.join(os.environ['CHARMM'],'toppar')
# Charmm executable:
charmm = os.path.join(os.environ['CHARMM'],'bin','charmm')
# construct the commandline
charmm_args = "HEADER=%(HEADER)s PSF=%(PSF)s DCD=%(DCD)s REF_PDB=%(REF_PDB)s "\
    "DCD_RMS=%(DCD_RMS)s "\
    "RECENTER=1" % job.filenames
inp = os.path.join(os.environ['CHARMM'],'analysis','trajectory','rmsfit.inp')
log = job.filenames['CHARMM_LOG']   # put name in locals()
cmd = "%(charmm)s %(charmm_args)s  <%(inp)s | tee %(log)s" % locals()
print cmd
os.system(cmd)


job.unstage()
job.cleanup()
