#!/usr/bin/env python
#$ -N rmsfit1IFCgsbp
#$ -S /usr/bin/python
#$ -v PYTHONPATH=/home/oliver/Library/python-lib
#$ -v LD_LIBRARY_PATH=/opt/intel/cmkl/8.0/lib/32:/opt/intel/itc60/slib:/opt/intel/ipp41/ia32_itanium/sharedlib:/opt/intel/ipp41/ia32_itanium/sharedlib/linux32:/opt/intel/fc/9.0/lib:/opt/intel/cc/9.0/lib
#$ -r n
#$ -j y
#$ -cwd
# $Id$

"""RMS-fit 1IFC GSBP trajectory to the I-FABP reference frame, based on a
STAMP alignment.

Note:

* This requires a centered trajectory in which the protein has been
  moved to the centre of the box. In Charmm do a

    'merge ...  recenter  select protein end ...'.

* Use the centered trajectory to fit on the reference frame.

However, it is much faster to do the following:

* Recenter (as above).
* Only fit the first frame of the recentered trajectory and write it out a crd.
* Use Charmm's 'merge ... orient' with the crd as reference. This is MUCH faster than
  using MDAnalysis for the rms fitting (20 Min vs 1-2 h).

"""

from staging.SunGridEngine import Job
#from staging.Local import Job
import os.path

#top_dir = '/nfs/animal/scratch0/oliver'
top_dir = '/Users/oliver/Biop/Projects'

seq_fasta = os.path.join(top_dir,'FABP','I-FABP','1IFC','GSBP','NVT','sequences','1IFC_GSBP15.fasta')

ref_dir = os.path.join(top_dir,'FABP','I-FABP','1IFC','coord')
ref_psf = os.path.join(ref_dir,'1ifc_xtal.psf')
ref_pdb = os.path.join(ref_dir,'1ifc_xtal.pdb')

trj_dir = os.path.join(top_dir,'FABP','I-FABP','1IFC','GSBP','NVT')
trj_psf = os.path.join(trj_dir,'inp','ifabp_gsbp_15_0.psf')
trj_dcd = os.path.join(trj_dir,'trj','ifabp_gsbp_15_1.dcd')

#------------------------------------------------------------
# EDIT THE inputfiles AND outputfiles DICTIONARIES.
#------------------------------------------------------------
# record input and output files relative to top_dir = cwd
job = Job(inputfiles=dict(sequence=seq_fasta,
                          ref_psf=ref_psf,
                          ref_pdb=ref_pdb,
                          trj_psf=trj_psf,
                          trj_dcd=trj_dcd,
                          ),
          outputfiles=dict(fit_dcd='trj/rmsfit_*.dcd',
                           ))
#
#------------------------------------------------------------

job.stage()

from MDAnalysis import Universe
import hop.trajectory

print "Alignment and selection string..."
selection = hop.trajectory.fasta2select(seq_fasta,is_aligned=True)

print "Setting up the Universes..."
ref = Universe(job.filenames['ref_psf'],pdbfilename=job.filenames['ref_pdb'])
trj = Universe(job.filenames['trj_psf'],job.filenames['trj_dcd'])

print "Fitting trajectory to reference..."
hop.trajectory.RMS_fit_trj(trj,ref, select=selection, prefix='rmsfit_')

job.unstage()
job.cleanup()
