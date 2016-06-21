#!/usr/bin/env python
#---------------- EDIT JOB NAME -------------------
#$ -N equivsites
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
job = Job(inputfiles=dict(dens_apo  = '1IFC/analysis/water.pickle',
                          dens_holo = '2IFB/analysis/water.pickle'),
          outputfiles=dict(saved_apo  = '1IFC/analysis/water.pickle',
                           saved_holo = '2IFB/analysis/water_remapped_1IFC.pickle',
                           graph = 'figs/equivalence_graph.png',
                           sitestats_apo = '1IFC/figs/sitestats.pdf',
                           sitestats_holo = '2IFB/figs/sitestats.pdf'))
#
#------------------------------------------------------------

job.stage()

# commands
import hop.utilities    # must come first (for reasons unknown)
hop.utilities.matplotlib_interactive(False)  # no X11 available

import hop.sitemap

class Struct:
    pass

dens = Struct()
dens.apo  = hop.sitemap.Density(filename=job.input['dens_apo'])    # reference
dens.holo0 = hop.sitemap.Density(filename=job.input['dens_holo'])  # will be remapped

print "=== remap with apo as reference ==="
dens.holo = hop.sitemap.remap_density(dens.holo0,dens.apo,verbosity=3)

print "=== find equivalence sites ==="
dens.holo.find_equivalence_sites_with(dens.apo,fmt='%d*',
                                      equivalence_graph=job.output['graph'],
                                      verbosity=7) # produces png at verb 7

print "=== ... and save the results ==="
dens.apo.save()
dens.holo.save(job.output['saved_holo'])

print "=== sitestats ==="
import hop.interactive
for density,fig in zip([dens.apo,dens.holo],['sitestats_apo','sitestats_holo']):
    figure = job.filenames[fig]
    hop.interactive.analyze_density(density,figure)


job.unstage()
job.cleanup()
