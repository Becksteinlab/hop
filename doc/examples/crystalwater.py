#!/usr/bin/env python
#---------------- EDIT JOB NAME -------------------
#$ -N xwat1OPA
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
""" Build a water density from crystal waters in the crystal structure."""

#from staging.SunGridEngine import Job
from staging.Local import Job

job = Job(inputfiles=dict(psf='',     # crystal psf
                          xtalpdb='', # crystal pdb (rms-fitted to reference)
                          refdensity='',  # reference water density
                          ),
          outputfiles=dict(xtal_density = 'analysis/xtal_water.pickle',
                           remapped_xtal = 'analysis/xtal_water_remapped_REFERENCE.pickle',
                           dx='analysis/*.dx')
          )
job.stage()
F = job.filenames
V = job.variables

import hop.utilities
hop.utilities.matplotlib_interactive(False)
import hop.sitemap

dens = hop.density.density_from_pdb(\
    psf=F['psf'], pdb=F['xtalpdb'],
    atomselection='segid XWAT and name OH2',
    delta=1.0,padding=4.0)

print "Mapping sites..."
dens.convert_density('water')
dens.map_sites(1.2)     # arbitrary cutoff...

print "Saving results..."
dens.save(F['xtal_density'])
dens.export()

print "Remapping to the reference density..."
ref = hop.sitemap.Density(filename=F['refdensity'])
remapped = hop.sitemap.remap_density(dens,ref,verbosity=3)
remapped.save(F['remapped_xtal'])
remapped.export()

job.unstage()
job.cleanup()
