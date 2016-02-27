#!/usr/bin/env python
#---------------- EDIT JOB NAME -------------------
#$ -N TAPdens
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
job = Job(inputfiles=dict(psf = 'inp/XXX.psf',
                          dcd = 'trj/rmsfit_XXX.dcd',
                          bulkdensity = 'analysis/bulk.pickle', # remove if needs computing
                          ),
          outputfiles=dict(TAP_dcd = 'trj/TAP_rmsfit_XXX.dcd',
                           dx = 'analysis/*.dx', 
                           waterdensity = 'analysis/water_TAP.pickle',
                           bulkdensity_saved = 'analysis/bulk.pickle'
                           ),
          variables = dict(water_selection = "name OH2",
                           rho_cut = 2.72,
                           rho_cut_bulk = 0.6,
                           )
          )
#
#------------------------------------------------------------

job.stage()
F = job.filenames
V = job.variables

# commands
import hop.utilities
hop.utilities.matplotlib_interactive(False)
from hop.interactive import *
import hop.sitemap, hop.trajectory
import MDAnalysis

# put first so that we fail soon if something is missing
print "== Bulk density =="
if 'bulkdensity' in job.input:
    print "Loading bulk density %s." % F['bulkdensity']
    bulk = hop.sitemap.Density(filename=F['bulkdensity'])
else:
    print "Computing bulk density %s." % F['bulkdensity_saved']
    print "Create bulk density from standard trajectory %s." % F['dcd']
    bulk = make_density(F['psf'],F['dcd'], F['bulkdensity_saved'],
                        atomselection=V['water_selection'],
                        soluteselection='protein and not name H*',cutoff=4)
bulk.map_sites(V['rho_cut_bulk'])


print "== TAP trajectory =="
u = MDAnalysis.Universe(F['psf'],F['dcd'])
woxy = u.select_atoms(V['water_selection'])
TAP = hop.trajectory.TAPtrajectory(trajectory=u.dcd,group=woxy,
                                   TAPradius=2.8,TAPsteps=3)
TAP.write(F['TAP_dcd'])

print "== Density =="
print "Create density from TAP trajectory %s." % F['TAP_dcd']
density = make_density(F['psf'],F['TAP_dcd'], F['waterdensity'],
                       atomselection=V['water_selection'])
density.map_sites(V['rho_cut'])
density.save()
density.export_map(combined=True)


print "== Combining densities =="
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

bulk.save()

density.save()
density.export_map(labels=1)  # write out dx file for bulk site

print "== Finishing... =="

job.unstage()
job.cleanup()
