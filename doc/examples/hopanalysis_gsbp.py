#!/usr/bin/env python
#---------------- EDIT JOB NAME -------------------
#$ -N hopanalysis
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
"""Complete analysis, starting from the density and the reference
density. This script is specialized for the GSBP setup (see below for
the usage using 'qsub -v GSBP_MD_DIR=...).

For convenience, use a bash script such as 

----------------------------------------------------------------------
#!/bin/bash
radius=15
ids="0 1 2"
function qsubscript () {
  local state=$1 id=$2
  qsub -v GSBP_MD_DIR=${state}/GSBP_MD_$id -N HA_R${radius}_${id}_${state} sge/hopanalysis_gsbp.py
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
# special setup for my stupid GSBP organization
import os,glob
import re
_statedir = {'apo':'1IFC', 'holo':'2IFB'}

# need to input full top dir in order to get (a) STATE, and (b) ID
try:
    GSBP_MD_DIR = os.environ['GSBP_MD_DIR']
except KeyError:
    raise KeyError('use -v GSBP_MD_DIR=1IFC/GSBP_MD_1 to specify the trajectory dir')
m = re.search('.*_MD_(?P<ID>[0-9]+).*',GSBP_MD_DIR)
try:
    ID = m.group('ID')    # ID is the number in the directory name (without _)
except:
    ID = ''
topdir = GSBP_MD_DIR.split(os.path.sep)[0]
_state = dict( [(dir,state) for state,dir in _statedir.items()] )
try:
    STATE = _state[topdir]  # calculate a apo or holo trajectory?
except KeyError:
    raise ValueError('STATE = '+str(STATE)+' must be one of '+str(_statedir.keys())+'.')

def TOP_DIR(*args):
    return os.path.join(topdir,*args)
def TRJ_DIR(*args):
    return os.path.join(GSBP_MD_DIR,*args)
def ANALYSIS_DIR(*args):
    return os.path.join(topdir,'analysis',*args)
def FIGS_DIR(*args):
    return os.path.join(topdir,'figs',*args)

def APO_DIR(*args):
    return os.path.join(_statedir['apo'],*args)
def HOLO_DIR(*args):
    return os.path.join(_statedir['holo'],*args)

psf_files = glob.glob(TRJ_DIR('*_gsbp_*.psf'))
psf_files.sort()
#ifabp_apo_gsbp_15_0.psf ifabp_apo_gsbp_gcmc_1.psf
psf = psf_files[-1]

dcd_files = glob.glob(TRJ_DIR('rmsfit_*_1.dcd'))  # note: _1 is NOT ID but refers to stage!
dcd_files.sort()
dcd = dcd_files[-1]

survival_times_dir = FIGS_DIR('survival_times_'+ID)

print "Input and Guesses for important variables:\n"\
    "\tSTATE=%(STATE)s  topdir=%(topdir)s\n"\
    "\tpsf=%(psf)s\n"\
    "\tdcd=%(dcd)s" % locals()

job = Job(variables=dict(trj=STATE,  # compute hop traj for trj = apo | holo HOP_TARGET
                         survival_times_dir=survival_times_dir,
                         ),
          inputfiles=dict(dens_apo  = APO_DIR('analysis','water_'+ID+'.pickle'),   # is modified!
                          dens_holo = HOLO_DIR('analysis','water_'+ID+'.pickle'),
                          psf = psf,
                          dcd = dcd,
                          ),
          outputfiles=dict(equiv_graph = FIGS_DIR('equivalence_graph_1.png'),
                           saved_apo = APO_DIR('analysis','*.pickle'),
                           saved_holo = HOLO_DIR('analysis','water_remapped_1IFC.pickle'),
                           hopdcd = TRJ_DIR('hoptrj.dcd'),
                           hoppsf = TRJ_DIR('hoptrj.psf'),
                           hopgraph = ANALYSIS_DIR('hopgraph_'+ID+'.pickle'),
                           xgmml = ANALYSIS_DIR('*.xgmml'),
                           survivaltimes=os.path.join(survival_times_dir,'*.png'),
                           ),
          )   
#
#------------------------------------------------------------
job.stage()
F = job.filenames
V = job.variables

# commands
import hop.utilities    # must come first (for reasons unknown)
hop.utilities.matplotlib_interactive(False)  # no X11 available

import hop.sitemap

dens = {}
dens['apo']  = hop.sitemap.Density(filename=F['dens_apo'])    # reference
dens['ref'] = dens['apo']
dens['holo0'] = hop.sitemap.Density(filename=F['dens_holo'])  # will be remapped

print "== Equivalence sites =="
print "=== remap with apo as reference ==="
if dens['holo0'].map.shape != dens['ref'].map.shape:
    print "Remapping for density "+str(dens['holo0'])+" required."
    dens['holo'] = hop.sitemap.remap_density(dens['holo0'],dens['ref'],verbosity=3)
else:
    print "Densities are already defined on the same grid, excellent."
    dens['holo'] = dens['holo0']

print "=== find equivalence sites ==="
dens['holo'].remove_equivalence_sites()   # start fresh
dens['ref'].remove_equivalence_sites()
dens['holo'].find_equivalence_sites_with(dens['ref'],fmt='%d*',
                                      equivalence_graph=F['equiv_graph'],
                                      verbosity=7) # produces png at verb 7

print "=== ... and save the results ==="
dens['apo'].save()                       # both modified densities ...
dens['holo'].save(F['saved_holo'])       # ... are needed later


print "== Hopping trajectory =="
import hop.interactive
import hop.constants

density = dens[V['trj']]  # select which density to make the hopping trajectory from

density.metadata['psf'] = F['psf'] # used in make_hoppingtraj()
density.metadata['dcd'] = F['dcd'] # (hidden API... bad practice :-p )

# Fixing metadata -- only necessary if the dcd has a wrong header
# as for instance produced by catdcd.
##delta_ps = 0.5   # the time between two frames in ps (Hirsh's trajectories)
##delta_AKMA = delta_ps * hop.constants.get_conversion_factor('time','ps','AKMA')
##density.metadata['delta'] = delta_AKMA
##fixtrajectory={'delta':density.metadata['delta']}
fixtrajectory = None

if fixtrajectory:
    print "Fixing new hopping trajectory with %r" % fixtrajectory

hoppingtrajectory = hop.interactive.make_hoppingtraj(density,F['hopdcd'],
                                                     fixtrajectory=fixtrajectory)


print "== Hopping graph =="

tgraph = hop.interactive.build_hoppinggraph(hoppingtrajectory,density)

h = tgraph.hopgraph       # main result is the 'hopgraph'
h.save(F['hopgraph'])     # save result
h.export(format='XGMML',use_filtered_graph=False) # as XGMML file, too

print "=== hop grap data ==="
h.filter(exclude={'outliers':True, 'Nmin':3, 'unconnected':True})
h.tabulate_k()            # show all calculated rate constants (filtered graph)
h.plot_fits(directory=V['survival_times_dir'])  # plot rate constant fits

print "== DONE =="
job.unstage()
job.cleanup()
