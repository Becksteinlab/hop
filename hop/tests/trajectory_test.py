# $Id$
# Specific profiling and testing code (I-FABP equilibrium simulation in
# /Users/oliver/Biop/Projects/WaterNetworks/testcases/1IFC )

def _setup_testcase_fromdensity():
    import MDAnalysis, hop.trajectory, hop.sitemap
    u = MDAnalysis.Universe('./NPT/ifabp_water.psf','./NPT/ifabp_water_1.dcd')
    group = u.select_atoms("name OH2")
    wd = hop.sitemap.Density(filename='waterdensity') 
    return hop.trajectory.HoppingTrajectory(u.dcd,group,wd)

def _setup_testcase_fromdcd():
    import hop.trajectory
    return hop.trajectory.HoppingTrajectory(hoppsf='hop.psf',hopdcd='hop.dcd')
            
def _test0(stop=2):
    """Check that the code is working"""
    sh = _setup_testcase_fromdensity()
    for ts in sh.map_dcd(0,stop):
        msg(3,"Mapping frame %5d/%d  [%5.1f%%]\r" % \
            (ts.frame,sh.n_frames,100.0*ts.frame/sh.n_frames))
        #print sh.frame
    msg(3,'\nMapping completed\n')
    return sh

def _test1(stop=5,repeat=5,N=4):
    """Timing map_dcd() --- see _make_timer()

    Canonical timing test when using the defaults

    16.2s   Id: sitehop.py 1135 2007-09-18 19:38:27Z oliver .. Heavy load?
     5.0s   $Id$: directly yield in map_dcd (gives 0.2 s), 11.4s under load
     
    """
    import timeit

    global _timer_mapdcd
    def _timer_mapdcd(sh,stop):
        for ts in sh.map_dcd(0,stop):
            msg(3,"Mapping frame %5d/%d  [%5.1f%%]\r" % \
                (ts.frame,sh.n_frames,100.0*ts.frame/sh.n_frames))
        msg(3,'\nMapping completed\n')


    tt = timeit.Timer(stmt='_timer_mapdcd(sh,'+str(stop)+')',
                      setup='from sitehop import _timer_mapdcd,_setup_testcase_fromdensity; '+\
                      'sh = _setup_testcase_fromdensity()')
    t = tt.repeat(N,repeat)
    tmin = min(t)
    print "All values: "+str(t)
    print "Min: "+str(tmin)+"s"
    return tmin

def _test2(stop=None,repeat=5,N=4):
    """Timing reading the trajectory

    Canonical timing test when using the defaults
    (Note: this reads the whole trajectorym not just 5 frames. stop is ignored)
    """
    import timeit

    global _timer_readhops
    def _timer_readhops(sh,stop):
        for ts in sh:
            msg(3,"Reading frame %5d/%d  [%5.1f%%]\r" % \
                (ts.frame,sh.n_frames,100.0*ts.frame/sh.n_frames))
        msg(3,'\nReading completed\n')
        
    tt = timeit.Timer(stmt='_timer_readhops(sh,'+str(stop)+')',
                      setup='from sitehop import _timer_readhops,_setup_testcase_fromdcd; '+\
                      'sh = _setup_testcase_fromdcd()')
    t = tt.repeat(N,repeat)
    tmin = min(t)
    print "All values: "+str(t)
    print "Min: "+str(tmin)+"s"
    return tmin
