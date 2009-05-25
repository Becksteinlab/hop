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
"""Complete analysis, starting from the density and the reference density.

Set the state for which trajectories and graphs are to be computed
from the qsub command line:

  qsub -v STATE=holo -N hopanalHOLO sge/hopanalysis.py

NOTE:
* currently set up to analyse I-FABP runs
* modify _statedir[] to customize for, eg CRBPII runs 
* will not work for GSBP runs (See hopanalysis_gsbp.py)
* enforces MINsite=1 (and rebuilds the densities; also saves remapped bulk; this
  bulk density is then defined on the same grid as the 'apo' or 'holo0' density,
  typically names 'analysis/water.pickle'

ASSUMPTIONS:

Bad heuristics to find psf and dcd:
* finds psf as last of 'inp/*.psf'
* finds dcd as last of 'trj/rmsfit_*-[0-9][0-9]ns.dcd'

Hard-coded paths for densities:
* apo and holo densities are under <APO/HOLO>/analysis/water.pickle
* apo and holo bulk densities are under <APO/HOLO>/analysis/bulk.pickle
"""
        
from staging.SunGridEngine import Job

#------------------------------------------------------------
# EDIT THE inputfiles AND outputfiles DICTIONARIES.
#------------------------------------------------------------
# record input and output files relative to top_dir = cwd
import os,glob
_statedir = {'apo':'1IFC', 'holo':'2IFB'}

try:
    STATE = os.environ['STATE']
except KeyError:
    raise KeyError('use -v STATE=apo to specify the state; script guesses the directory')

try:
    topdir = _statedir[STATE]  # calculate a apo or holo trajectory?
except KeyError:    
    raise ValueError('STATE = '+str(STATE)+' must be one of '+str(_statedir.keys())+'.')

def TOP_DIR(*args):
    return os.path.join(topdir,*args)
def TRJ_DIR(*args):
    return os.path.join(topdir,'trj',*args)
def ANALYSIS_DIR(*args):
    return os.path.join(topdir,'analysis',*args)
def FIGS_DIR(*args):
    return os.path.join(topdir,'figs',*args)

def APO_DIR(*args):
    return os.path.join(_statedir['apo'],*args)
def HOLO_DIR(*args):
    return os.path.join(_statedir['holo'],*args)

psf_files = glob.glob(TOP_DIR('inp','*.psf'))                # CHECK THIS
psf_files.sort()
# inp/ifabp_apo_100mM_charmm_atoms.psf  inp/ifabp_apo_100mM_charmm.psf
psf = psf_files[-1]

dcd_files = glob.glob(TRJ_DIR('rmsfit_*-[0-9][0-9]ns.dcd'))  # CHECK THIS
dcd_files.sort()
# trj/rmsfit_ifabp_apo_100mM_0-1ns.dcd  trj/rmsfit_ifabp_apo_100mM_0-20ns.dcd
dcd = dcd_files[-1]

survival_times_dir = FIGS_DIR('survival_times')

print "Input and Guesses for important variables:\n"\
    "\tSTATE=%(STATE)s  topdir=%(topdir)s\n"\
    "\tpsf=%(psf)s\n"\
    "\tdcd=%(dcd)s" % locals()

job = Job(variables=dict(trj=STATE,  # compute hop traj for trj = apo | holo HOP_TARGET
                         survival_times_dir=survival_times_dir,
                         MINsite=1,
                         ),
          inputfiles=dict(dens_apo  = APO_DIR('analysis','water.pickle'),   # is modified!
                          dens_holo = HOLO_DIR('analysis','water.pickle'),
                          bulk_apo = APO_DIR('analysis','bulk.pickle'),
                          bulk_holo = HOLO_DIR('analysis','bulk.pickle'),
                          psf = psf,
                          dcd = dcd,
                          ),
          outputfiles=dict(equiv_graph = FIGS_DIR('equivalence_graph.png'),
                           saved_apo = APO_DIR('analysis','*.pickle'),
                           saved_holo = HOLO_DIR('analysis','water_remapped_1IFC.pickle'),
                           hopdcd = TRJ_DIR('hoptrj.dcd'),
                           hoppsf = TRJ_DIR('hoptrj.psf'),
                           hopgraph = ANALYSIS_DIR('hopgraph.pickle'),
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

print "== Loading densities =="
dens = {}
dens['apo']   = hop.sitemap.Density(filename=F['dens_apo'])    # reference
dens['ref']   = dens['apo']
dens['holo0'] = hop.sitemap.Density(filename=F['dens_holo'])  # will be remapped
bulkdens = {}
bulkdens['apo']   = hop.sitemap.Density(filename=F['bulk_apo'])
bulkdens['holo0']  = hop.sitemap.Density(filename=F['bulk_holo'])

print "Loaded densities %r" % dens.items()
print "Loaded bulk densities %r" % bulkdens.items()

def map_sites(state,rhocut=2.72):
    _prefx = "[%(state)s] density: " % locals()
    density = dens[state]
    bulk    = bulkdens[state]    
    density.P['MINsite'] = V['MINsite']    # necessary for older densities
    density.map_sites(rhocut)              # erases bulk

    print _prefx + "mapped sites for MINsite=%(MINsite)d at rho_cut=%(threshold)f" % density.P
    print _prefx + "adding the biggest bulk site at position 1" % locals()
    try:
        density.site_insert_bulk(bulk)   # hack!
    except ValueError,msg:
        print msg
        print "Remapping bulk onto density..."
        import hop.sitemap
        bulk_remapped = hop.sitemap.remap_density(bulk,density,verbosity=3)
        bulk_remapped.filename(bulk.filename(),set_default=True)  # we will overwrite original bulk        
        bulk = bulk_remapped
        density.site_insert_bulk(bulk)   # hack!
    density.save()
    bulk.save()

print "Recalculating the site map with MINsite=%(MINsite)d" % V
map_sites('apo')
map_sites('holo0')
del bulkdens

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
delta_ps = 1.0   # the time between two frames in ps (LAMMPS)
delta_AKMA = delta_ps * hop.constants.get_conversion_factor('time','ps','AKMA')
density.metadata['delta'] = delta_AKMA
fixtrajectory={'delta':density.metadata['delta']}
# fixtrajectory = None

if fixtrajectory:
    print "Fixing new hopping trajectory with %r" % fixtrajectory

hoppingtrajectory = hop.interactive.make_hoppingtraj(density,F['hopdcd'],
                                                     fixtrajectory=fixtrajectory)


print "== Hopping graph =="

tgraph = hop.interactive.build_hoppinggraph(hoppingtrajectory,density)

h = tgraph.hopgraph       # main result is the 'hopgraph'
h.save(F['hopgraph'])     # save result
h.export(format='XGMML',use_filtered_graph=False) # as XGMML file, too

print "=== hop graph data ==="
h.filter(exclude={'outliers':True, 'Nmin':3, 'unconnected':True})
h.tabulate_k()            # show all calculated rate constants (filtered graph)
h.plot_fits(directory=V['survival_times_dir'])  # plot rate constant fits

print "== DONE =="
job.unstage()
job.cleanup()
