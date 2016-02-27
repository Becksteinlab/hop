#!/usr/bin/env python
#---------------- EDIT JOB NAME -------------------
#$ -N SiteAnalysis
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
"""Analyse a trajectory in terms of properties of hydration sites.

This script is specialized for the GSBP setup (see below for the usage
using 'qsub -v GSBP_MD_DIR=...).

For convenience, use a bash script:

----------------------------------------------------------------------
#!/bin/bash
radius=15
ids="0 1 2"
function qsubscript () {
  local state=$1 id=$2
  qsub -v GSBP_MD_DIR=${state}/GSBP_MD_$id -N SA_R${radius}_${id}_${state} sge/siteanalysis_gsbp.py
}

for state in 1IFC 2IFB; do
  for id in $ids; do
    echo ">> state=$state id=$id"
    qsubscript $state $id
  done
done
----------------------------------------------------------------------
"""
        
from staging.SunGridEngine import Job

#------------------------------------------------------------
# EDIT THE inputfiles AND outputfiles DICTIONARIES.
#------------------------------------------------------------
# record input and output files relative to top_dir = cwd
import os
import re,glob
try:
    GSBP_MD_DIR = os.environ['GSBP_MD_DIR']
except KeyError:
    raise KeyError('use -v GSBP_MD_DIR=1IFC/GSBP_MD_1 to specify the trajectory dir')
m = re.search('.*_MD_(?P<id>[0-9]+).*',GSBP_MD_DIR)
try:
    id = m.group('id')    # id is the number in the directory name (without _)
except:
    id = ''
topdir = GSBP_MD_DIR.split(os.path.sep)[0]
if topdir not in ('1IFC','2IFB'):
    raise ValueError('topdir = '+str(topdir)+' must be one of 1IFC or 2IFB')

def TRJ_DIR(*args):
    return os.path.join(GSBP_MD_DIR,*args)

def TOP_DIR(*args):
    return os.path.join(topdir,*args)


psf_files = glob.glob(TRJ_DIR('*_gsbp_*.psf'))
psf_files.sort()
#ifabp_apo_gsbp_15_0.psf ifabp_apo_gsbp_gcmc_1.psf
psf = psf_files[-1]

dcd_files = glob.glob(TRJ_DIR('rmsfit_*_1.dcd'))
dcd_files.sort()
dcd = dcd_files[-1]

density_files = glob.glob(TOP_DIR('analysis','water_remapped*_'+id+'.pickle'))
if len(density_files) == 0:
    density_files = glob.glob(TOP_DIR('analysis','water_'+id+'.pickle'))
density_files.sort()
# 2IFB/analysis/water_1.pickle 2IFB/analysis/water_remapped_1IFC_1.pickle
density = density_files[-1]  # pick remapped if available

job = Job(variables=dict(),
          inputfiles=dict(psf=psf,
                          trj=dcd,
                          hoppsf=TRJ_DIR('hoptrj.psf'),
                          hoptrj=TRJ_DIR('hoptrj.dcd'),
                          density=density,
                          ),
          outputfiles=dict(fig_distance=TOP_DIR('figs','site_distance_'+id+'.pdf'),
                           fig_orbit=TOP_DIR('figs','site_orbit_'+id+'.pdf'),
                           fig_distance2=TOP_DIR('figs','site_distance2_'+id+'.pdf'),
                           fig_orbit2=TOP_DIR('figs','site_orbit2_'+id+'.pdf'),
                           fig_occupancy=TOP_DIR('figs','site_occupancy_'+id+'.pdf'),
                           fig_orbitoccupancy=TOP_DIR('figs','site_orbitoccupancy_'+id+'.pdf'),
                           saved=TOP_DIR('analysis','siteanalysis_'+id+'.pickle'),
                           ))
#
#------------------------------------------------------------
job.stage()
F = job.filenames
V = job.variables

# commands
import hop.utilities    # must come first (for reasons unknown)
hop.utilities.matplotlib_interactive(False)  # no X11 available

import MDAnalysis
import hop.siteanalysis, hop.sitemap, hop.trajectory

print "== loading data =="
u = MDAnalysis.Universe(F['psf'],F['trj'])
water = u.select_atoms('name OH2')
hoptraj = hop.trajectory.HoppingTrajectory(hoppsf=F['hoppsf'],hopdcd=F['hoptrj'])
density = hop.sitemap.Density(filename=F['density'])

A = hop.siteanalysis.Collection(water,hoptraj,density)

print "== computing observables for all sites =="
A.add_observable('distance',lo_bin=0,hi_bin=4)
A.add_observable('orbit',lo_bin=0,hi_bin=10)
A.add_observable('occupancy',lo_bin=0,hi_bin=4)
A.add_observable('orbitoccupancy',lo_bin=0,hi_bin=8)
A.compute()
A.save(F['saved'])
A.plot('distance',filename=F['fig_distance'])
A.plot('orbit',filename=F['fig_orbit'])
A.imshow('distance',vmax=None,filename=F['fig_distance2'])
A.imshow('orbit',vmax=None,filename=F['fig_orbit2'])
A.imshow('occupancy',filename=F['fig_occupancy'])
A.imshow('orbitoccupancy',filename=F['fig_orbitoccupancy'])

print "== DONE =="

job.unstage()
job.cleanup()
