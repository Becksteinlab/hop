#!/usr/bin/env python
#---------------- EDIT JOB NAME -------------------
#$ -N combgraph
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
# from staging.Local import Job

#------------------------------------------------------------
# EDIT THE inputfiles AND outputfiles DICTIONARIES.
#------------------------------------------------------------
# record input and output files relative to top_dir = cwd
job = Job(inputfiles=dict(hg_apo =  '1IFC/analysis/hopgraph.pickle',
                          hg_holo = '2IFB/analysis/hopgraph.pickle'),
          outputfiles=dict(cg = 'analysis/cg_1IFC_2IFB.pickle',
                           png = 'fig/*.png') )
#
#------------------------------------------------------------

job.stage()

# commands
import hop.utilities    # must come first (for reasons unknown)
hop.utilities.matplotlib_interactive(False)  # no X11 available

import hop.sitemap,hop.graph


class Struct:
    pass

hg = Struct()
hg.apo  = hop.graph.HoppingGraph(filename=job.input['hg_apo'])
hg.holo = hop.graph.HoppingGraph(filename=job.input['hg_holo'])

cg = hop.graph.CombinedGraph(g0=hg.apo,g1=hg.holo)
cg.filter(exclude={'outliers':True,'bulk':True, 'unconnected':False})
cg.save(job.output['cg'])

cg.plot(0,'fig/cg_apo',linewidths=(0.01,))   # produces png
cg.plot(1,'fig/cg_holo',linewidths=(0.01,))

job.unstage()
job.cleanup()
