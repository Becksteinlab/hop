#!/usr/bin/env python
#---------------- EDIT JOB NAME -------------------
#$ -N graph
#--------------------------------------------------
#$ -S /usr/bin/python
#$ -v PATH
#$ -v PYTHONPATH=/home/oliver/Library/python-lib
#$ -v LD_LIBRARY_PATH=/opt/intel/cmkl/8.0/lib/32:/opt/intel/itc60/slib:/opt/intel/ipp41/ia32_itanium/sharedlib:/opt/intel/ipp41/ia32_itanium/sharedlib/linux32:/opt/intel/fc/9.0/lib:/opt/intel/cc/9.0/lib
#$ -r n
#$ -j y
# Using the current working directory is IMPORTANT with the default settings for Job()
#$ -cwd
#$ -m e
# $Id$

"""Analyse a hopping trajectory in terms of a graph.

* creates hopping graph
* exports xgml file
* initial analysis and statistics
"""


from staging.SunGridEngine import Job

#------------------------------------------------------------
# EDIT THE inputfiles AND outputfiles DICTIONARIES.
#------------------------------------------------------------
# record input and output files relative to top_dir = cwd
job = Job(inputfiles=dict(hopdcd='trj/hops.dcd',
                          hoppsf='trj/hops.psf',
                          density='analysis/water.pickle',
                          ),
          outputfiles=dict(xgmmlfiles='analysis/*.xgmml',
                           hopgraph='analysis/hopgraph.pickle',
                           survivaltimes='analysis/survival_times/*.png',
                          ))
#
#------------------------------------------------------------

job.stage()

# commands
import hop.utilities
hop.utilities.matplotlib_interactive(False)  # no X11 available

from hop.interactive import build_hoppinggraph_fromfiles

tgraph = build_hoppinggraph_fromfiles(job.input['hoppsf'],job.input['density'])

# analysis of hop graph   
h = tgraph.hopgraph       # main result is the 'hopgraph'
h.save(job.output['hopgraph'])

print "== Tabulated rate constants =="
h.filter(exclude={'outliers':True, 'Nmin':3, 'unconnected':True})
h.tabulate_k()            # show all calculated rate constants (filtered graph)

print "== Detailed statistics for all sites =="
for site in h.graph.nodes():
    h.show_site(site,use_filtered_graph=False)
    print

print "== Plot fits =="
h.plot_fits(directory='analysis/survival_times')  # plot rate constant fits

print "== Export graphs as xgmml =="
print "* write xgmml file to visualize full graph"
h.export('analysis/hg_water',use_filtered_graph=False)

print "* write xgmml file to visualize (filtered) graph without bulk"
h.filter(exclude={'outliers':True,'bulk':True})
h.export('analysis/hg_water_nobulk')


job.unstage()
job.cleanup()   
