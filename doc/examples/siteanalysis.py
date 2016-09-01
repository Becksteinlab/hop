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
"""Analyse a trajectory in terms of properties of hydration sites."""
        
from staging.SunGridEngine import Job

#------------------------------------------------------------
# EDIT THE inputfiles AND outputfiles DICTIONARIES.
#------------------------------------------------------------
# record input and output files relative to top_dir = cwd
job = Job(variables=dict(),
          inputfiles=dict(psf="",
                          trj="",
                          hoppsf="",
                          hoptrj="",
                          density="",
                          ),
          outputfiles=dict(fig_distance="figs/site_distance.pdf",
                           fig_orbit="figs/site_orbit.pdf",
                           fig_distance2="figs/site_distance2.pdf",
                           fig_orbit2="figs/site_orbit2.pdf",
                           fig_occupancy="figs/site_occupancy.pdf",
                           fig_orbitoccupancy="figs/site_orbitoccupancy.pdf",
                           saved="analysis/siteanalysis.pickle",
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
A.add_observable('occupancy',lo_bin=0,hi_bin=4,trajectory=True)
A.add_observable('orbitoccupancy',lo_bin=0,hi_bin=10,trajectory=True)
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
