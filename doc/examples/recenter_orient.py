#!/usr/bin/env python
#---------------- EDIT JOB NAME -------------------
#$ -N NAME
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
#
# see http://gonzo.med.jhmi.edu/woolfwiki/index.php/Making_proteins_whole_in_LAMMPS_trajectories
#
        
from staging.SunGridEngine import Job

#------------------------------------------------------------
# EDIT THE inputfiles AND outputfiles DICTIONARIES.
#------------------------------------------------------------
# record input and output files relative to top_dir = cwd
job = Job(inputfiles=dict(HEADER='str/header.str',
                          PSF='inp/ifabp_apo_100mM_charmm.psf',
                          DCD='trj/ifabp_apo_100mM_0-20ns_dt1ns.dcd',    
                          REF_PDB='coord/1ifc_xtal.pdb',
                          ),
          outputfiles=dict(DCD_RMS='trj/rmsfit_ifabp_apo_100mM_0-20ns_dt1ns.dcd',
                           RESINDEX='inp/ifabp_apo_100mM_charmm.resindex',
                           LOG_1='log/1_remap.log',
                           LOG_2='log/2_unfold.log',
                           LOG_3='log/3_remap.log',
                           LOG_4='log/4_orient.log',
                           ))
#
#------------------------------------------------------------

job.stage()

#------------------------------------------------------------
# python script proper starts here
# (build commandlines using python named string interpolation)
#
import os

# This setup requires my canonical filesystem layout under CHARMM:
if not 'CHARMM' in os.environ:   # get sane value for CHARMM top dir
    os.environ['CHARMM'] = os.path.join(os.environ['BIO_D'],'Library','Charmm')
    os.environ['CHARMM_TOP'] = os.path.join(os.environ['CHARMM'],'toppar')
if not 'LAMMPS' in os.environ:   # get sane value for LAMMPS top dir
    os.environ['LAMMPS'] = os.path.join(os.environ['BIO_D'],'Library','LAMMPS')

# executables
executables = {'charmm':os.path.join(os.environ['CHARMM'],'bin','charmm'),
               'remap_dcd':os.path.join(os.environ['LAMMPS'],'bin','trjconv_dcd'),
               'extract_resindex':os.path.join(os.environ['LAMMPS'],'bin','extract_resindex.py'),
               }

# dict of files to build command lines from using named interpolation
inp = {'rmsfit':os.path.join(os.environ['CHARMM'],'analysis','trajectory','rmsfit.inp'),
       'DCD_1':'trj/1.dcd',
       'DCD_2':'trj/2.dcd',
       'DCD_3':'trj/3.dcd',
       'IMAGE_CUTOFF':50,
       'RES_1':1, 'RES_2':131, # define group to recenter on: RES_1...RES_2 (index numbers!)
       }
inp.update(executables)
inp.update(job.filenames)

# 0. build resindex
cmd = "%(extract_resindex)s %(PSF)s %(RESINDEX)s" % inp
print cmd
os.system(cmd)

# 1. remap waters
cmd = "%(remap_dcd)s --index=%(RESINDEX)s --outdcd=%(DCD_1)s %(DCD)s  >%(LOG_1)s" % inp
print cmd
os.system(cmd)

# 2. unfold and recenter
cmd = "%(charmm)s <%(rmsfit)s  >%(LOG_2)s "\
    "HEADER=%(HEADER)s PSF=%(PSF)s DCD=%(DCD_1)s "\
    "DCD_RMS=%(DCD_2)s "\
    "IMAGE_CUTOFF=%(IMAGE_CUTOFF)d "\
    "RECENTER=1 UNFOLD=1 ORIENT=0" % inp
print cmd
os.system(cmd)
os.unlink(inp['DCD_1'])  # clean up space

# 3. remap water again
cmd = "%(remap_dcd)s --index=%(RESINDEX)s --outdcd=%(DCD_3)s "\
    "--pbc=atoms --center=resids --r1=%(RES_1)d --r2=%(RES_2)d "\
    "%(DCD_2)s >%(LOG_3)s" % inp
print cmd
os.system(cmd)
os.unlink(inp['DCD_2'])  # clean up space

# 4. orient (RMS-fit) and recenter a few stray waters
cmd = "%(charmm)s <%(rmsfit)s  >%(LOG_4)s "\
    "HEADER=%(HEADER)s PSF=%(PSF)s DCD=%(DCD_3)s "\
    "DCD_RMS=%(DCD_RMS)s "\
    "REF_PDB=%(REF_PDB)s "\
    "IMAGE_CUTOFF=50 "\
    "RECENTER=1 UNFOLD=0 ORIENT=1" % inp
print cmd
os.system(cmd)
os.unlink(inp['DCD_3'])  # clean up space


job.unstage()
job.cleanup()
