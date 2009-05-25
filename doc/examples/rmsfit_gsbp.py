#!/usr/bin/env python
#---------------- EDIT JOB NAME -------------------
#$ -N gsbpRMSfit
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
"""RMS-fit a GSBP trajectory (with fixed atoms and gaps in the
protein) to a reference structure.  This requires a sequence
alignment, the reference, and the trajectory.  All intermediate steps
such as creating the initial aligned frame (using MDAnalysis) are
performed, as is the actual fitting using Charmm with the CONSFIX=1
ORIENT=1 options.

Currently, this is geared towards GSBP_MD simulations. Set GSBP_MD_DIR in the environment or with

  qsub -v GSBP_MD_DIR=./GSBP_MD_1 rmsfit_gsbp.py

It automatically find the psf file for the dcd in GSBP_MD_DIR, but other values have to be set in 
Job() as usual.
"""
        
from staging.SunGridEngine import Job

#------------------------------------------------------------
# EDIT THE inputfiles AND outputfiles DICTIONARIES.
#------------------------------------------------------------
# record input and output files relative to top_dir = cwd
import os,glob
try:
    GSBP_MD_DIR = os.environ['GSBP_MD_DIR']
except KeyError:
    #raise KeyError('use -v GSBP_MD_DIR=./GSBP_MD_1 to specify the trajectory dir')
    GSBP_MD_DIR = './GSBP_MD_0'
def TRJ_DIR(*args):
    return os.path.join(GSBP_MD_DIR,*args)

psf_files = glob.glob(TRJ_DIR('*_gsbp_*.psf'))
psf_files.sort()
#ifabp_apo_gsbp_15_0.psf ifabp_apo_gsbp_gcmc_1.psf
psf = psf_files[-1]

job = Job(inputfiles=dict(HEADER=TRJ_DIR('header.str'),
                          PSF=psf,
                          DCD=TRJ_DIR('ifabp_apo_gsbp_15_1.dcd'),    
                          sequence='sequences/1IFC_R151IFC.fasta',
                          ref_psf='coord/1ifc_xtal.psf',
                          ref_pdb='coord/1ifc_xtal.pdb',
                          trj_psf='GSBPsetup/ifabp_apo_gsbp_15_0.psf',
                          trj_pdb='GSBPsetup/ifabp_apo_gsbp_15_0.pdb',
                          ),
          outputfiles=dict(REF_PDB='GSBPsetup/rmsfit_ifabp_apo_gsbp_15_0.pdb',
                           DCD_RMS=TRJ_DIR('rmsfit_ifabp_apo_gsbp_15_1.dcd'),
                           LOG_1='log/1_referenceframe.log',
                           LOG_2='log/2_rmsfit_charmm.log',
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
       }
inp.update(executables)
inp.update(job.filenames)


# 1. build reference frame from alignment
#------------------------------------------------------------
from MDAnalysis import Universe
import hop.trajectory

print "Setting up the Universes..."
ref = Universe(job.filenames['ref_psf'],pdbfilename=job.filenames['ref_pdb'])
trj = Universe(job.filenames['trj_psf'],job.filenames['trj_pdb'])
ref_resids    = [a.resid for a in ref.selectAtoms('name CA')]
target_resids = [a.resid for a in trj.selectAtoms('name CA')]
print "Alignment and selection string..."
selection = hop.trajectory.fasta2select(job.filenames['sequence'],
                                        ref_resids=ref_resids,target_resids=target_resids,
                                        is_aligned=True)
print "Fitting trajectory to reference..."
hop.trajectory.RMS_fit_trj(trj,ref, select=selection, filename=job.filenames['REF_PDB'])
print "Done: result is '%(REF_PDB)s'" % job.filenames
#------------------------------------------------------------

# 2. orient (RMS-fit)
print "RMS-fitting using Charmm"
cmd = "%(charmm)s <%(rmsfit)s  >%(LOG_2)s "\
    """HEADER='"%(HEADER)s"' PSF='"%(PSF)s"' DCD='"%(DCD)s"' """\
    """DCD_RMS='"%(DCD_RMS)s"' """\
    """REF_PDB='"%(REF_PDB)s"' """\
    """RECENTER=0 UNFOLD=0 ORIENT=1 CONSFIX=1""" % inp
print cmd
os.system(cmd)

job.unstage()
job.cleanup()
