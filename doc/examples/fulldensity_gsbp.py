#!/usr/bin/env python
#---------------- EDIT JOB NAME -------------------
#$ -N den1IFCg0
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
"""Make the density pickle files from a rms fitted trajectory.
Uses flexible GSBP_MD_DIR.

usage: qsub -v GSBP_MD_DIR=./GSBP_MD_1  fulldensity_gsbp.py

to specify the trajectory dir

"""
        
from staging.SunGridEngine import Job

#------------------------------------------------------------
# EDIT THE inputfiles AND outputfiles DICTIONARIES.
#------------------------------------------------------------
# record input and output files relative to top_dir = cwd
import os,glob
import re
try:
    GSBP_MD_DIR = os.environ['GSBP_MD_DIR']
except KeyError:
    #raise KeyError('use -v GSBP_MD_DIR=./GSBP_MD_1 to specify the trajectory dir')
    GSBP_MD_DIR = './GSBP_MD_0'
m = re.search('.*_MD_(?P<id>[0-9]+).*',GSBP_MD_DIR)
try:
    id = m.group('id')    # id is the number in the directory name (without _)
except:
    id = ''
def TRJ_DIR(*args):
    return os.path.join(GSBP_MD_DIR,*args)

psf_files = glob.glob(TRJ_DIR('*_gsbp_*.psf'))
psf_files.sort()
#ifabp_apo_gsbp_15_0.psf ifabp_apo_gsbp_gcmc_1.psf
psf = psf_files[-1]

dcd_files = glob.glob(TRJ_DIR('rmsfit_*_1.dcd'))
dcd_files.sort()
dcd = dcd_files[-1]

job = Job(inputfiles=dict(psf = psf,
                          dcd = dcd),
          outputfiles=dict(dx = 'analysis/*.dx', 
                           water_saved = 'analysis/water_'+id+'.pickle',
                           bulk_saved = 'analysis/bulk_'+id+'.pickle'))
#
#------------------------------------------------------------

job.stage()

# commands
import hop.utilities
hop.utilities.matplotlib_interactive(False)
from hop.interactive import *

density = make_density(job.filenames['psf'],job.filenames['dcd'],job.filenames['water_saved'],
                       atomselection="name OH2")
density.map_sites(2.72)
density.save()
density.export_map(combined=True)

bulk = make_density(job.filenames['psf'],job.filenames['dcd'],job.filenames['bulk_saved'],
                    atomselection="name OH2",soluteselection='protein and not name H*',cutoff=4)
# VMD only works on small trajectories !!
#bulk = make_density(job.filenames['psf'],job.filenames['dcd'], './analysis/bulk',
#                    atomselection="name OH2 and not within 4 of (protein and not hydrogen)",
#                    backend='VMD')
bulk.map_sites(0.6)

# add the biggest bulk site at position 1
try:
    density.site_insert_bulk(bulk)   # hack!
except ValueError,msg:
    print msg
    print "Remapping bulk onto density..."
    import hop.sitemap
    bulk_remapped = hop.sitemap.remap_density(bulk,density,verbosity=3)
    bulk = bulk_remapped
    density.site_insert_bulk(bulk)   # hack!

bulk.save(job.filenames['bulk_saved'])

density.save()
density.export_map(labels=1)  # write out dx file for bulk site

job.unstage()
job.cleanup()
