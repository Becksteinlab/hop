"""Align coordinate to reference, using a fasta alignment to find the
common atoms."""

from staging.SunGridEngine import Job
#from staging.Local import Job

job = Job(inputfiles=dict(sequence='sequences/1IFC_R151IFC.fasta',
                          ref_psf='coord/1ifc_xtal.psf',
                          ref_pdb='coord/1ifc_xtal.pdb',
                          trj_psf='GSBPsetup/ifabp_apo_gsbp_15_0.psf',
                          trj_pdb='GSBPsetup/ifabp_apo_gsbp_15_0.pdb',
                          ),
          outputfiles=dict(fit_pdb='GSBPsetup/rmsfit_ifabp_apo_gsbp_15_0.pdb',
                           ))

job.stage()

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
hop.trajectory.RMS_fit_trj(trj,ref, select=selection, filename=job.filenames['fit_pdb'])

print "Done: result is '%(fit_pdb)s'" % job.filenames

job.unstage()
job.cleanup()
