#!/usr/bin/env python
#---------------- EDIT JOB NAME -------------------
#$ -N density
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
job = Job(inputfiles=dict(psf = '1IFC/inp/XXX.psf',
                          dcd = '1IFC/trj/rmsfit_XXX.dcd'),
          outputfiles=dict(water = '1IFC/analysis/water.pickle',
                           bulk = '1IFC/analysis/bulk.pickle',
                           dx = '1IFC/analysis/*.dx',
                           ))
#
#------------------------------------------------------------

job.stage()

F = job.filenames
V = job.variables

# commands
import hop.utilities
hop.utilities.matplotlib_interactive(False)
from hop.interactive import *

hop.utilities.set_verbosity(3)

print "== making density =="
density = make_density(F['psf'],F['dcd'],F['water'],atomselection="name OH2")
density.map_sites(2.72)
density.save()
density.export_map(combined=True)

print "== making bulk density =="
bulk = make_density(F['psf'],F['dcd'], F['bulk'],
         atomselection="name OH2",soluteselection='protein and not name H*',cutoff=4)
# VMD only works on small trajectories !!
#bulk = make_density(job.input['psf'],job.input['dcd'], './analysis/bulk',
#                    atomselection="name OH2 and not within 4 of (protein and not hydrogen)",
#                    backend='VMD')
bulk.map_sites(0.6)

print "adding the biggest bulk site at position 1"
try:
    density.site_insert_bulk(bulk)   # hack!
except ValueError,msg:
    print msg
    print "Remapping bulk onto density..."
    import hop.sitemap
    bulk_remapped = hop.sitemap.remap_density(bulk,density,verbosity=3)
    bulk = bulk_remapped
    density.site_insert_bulk(bulk)   # hack!

print "== Saving results =="
bulk.save()
density.save()

#density.export_map(labels=1)  # write out dx file for bulk site

job.unstage()
job.cleanup()
