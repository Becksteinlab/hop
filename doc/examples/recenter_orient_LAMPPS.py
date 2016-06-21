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
# $Id: recenter_orient.py 4090 2010-04-08 19:38:34Z oliver $
#
# see http://gonzo.med.jhmi.edu/woolfwiki/index.php/Making_proteins_whole_in_LAMMPS_trajectories
#
"""Prepare LAMPPS trajectories that are suitable as input for density calculations:
1. centred on the protein
2. rms-fitted to reference
3. protein and water molecules are whole

This uses the new ``trjconv_dcd`` (~rev 120).
"""

#from staging.SunGridEngine import Job
from staging.Local import Job

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
                           LOG_2='log/2_orient.log',
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
               'trjconv_dcd':os.path.join(os.environ['LAMMPS'],'bin','trjconv_dcd'),
               'extract_resindex':os.path.join(os.environ['LAMMPS'],'bin','trjconv_resindex.py'),
               }

# dict of files to build command lines from using named interpolation
inp = {'rmsfit':os.path.join(os.environ['CHARMM'],'analysis','trajectory','rmsfit.inp'),
       'DCD_1':'trj/1.dcd',
       'IMAGE_CUTOFF':50,
       'RES_1':1, 'RES_2':131, # define group to recenter on: RES_1...RES_2 (index numbers!)
       }
inp.update(executables)
inp.update(job.filenames)

# 0. build resindex
cmd = "%(extract_resindex)s %(PSF)s %(RESINDEX)s > %(LOG_0)s" % inp
print cmd
rc = os.system(cmd)
if rc != 0:
    raise OSError("System command %(cmd)r failed with error code %(rc)r" % vars())

# 1. remap waters and make protein whole
if not os.path.exists(inp['DCD_1']):
    cmd = """%(trjconv_dcd)s --index=%(RESINDEX)s \
                         --center=resids --r1=%(RES_1)d --r2=%(RES_2)d --remap_selection \
                         --pbc=atoms \
                         --outdcd=%(DCD_1)s %(DCD)s  >%(LOG_1)s""" % inp
    print cmd
    rc = os.system(cmd)
    if rc != 0:
        raise OSError("System command %(cmd)r failed with error code %(rc)r" % vars())
else:
    print "%(DCD_1)s already exists, skipping trjconv_dcd step" % inp

# 2. orient (RMS-fit) 
# XXX RECENTER=1: and recenter a few stray waters
cmd = """%(charmm)s <%(rmsfit)s  >%(LOG_2)s """\
    """HEADER='"%(HEADER)s"' PSF='"%(PSF)s"' DCD='"%(DCD_1)s"' """\
    """DCD_RMS='"%(DCD_RMS)s"' """\
    """REF_PDB='"%(REF_PDB)s"' """\
    'IMAGE_CUTOFF=%(IMAGE_CUTOFF)f '\
    'RECENTER=0 UNFOLD=0 ORIENT=1' % inp
print cmd
rc = os.system(cmd)
if rc != 0:
    raise OSError("System command %(cmd)r failed with error code %(rc)r" % vars())

os.unlink(inp['DCD_1'])  # clean up space


job.unstage()
job.cleanup()
