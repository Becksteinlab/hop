import hop.sitemap
import MDAnalysis
import numpy as np
import copy

def flux_graph(topology, density,):
    
    u=MDAnalysis.Universe(topology)
    bins=[l for l in xrange(u.dimensions[2])]
    p=u.selectAtoms("protein")
    l=u.selectAtoms("resname POP")
    pmin,pmax=p.atoms.bbox()
    lmin,lmax=l.atoms.bbox()
    bmin,bmax=[0,0,0],[u.dimensions[0],u.dimensions[1],u.dimensions[2]]
    zmin_p,zmax_p=np.digitize([pmin[2],pmax[2]],bins)
    ymin_p,ymax_p=np.digitize([pmin[1],pmax[1]],bins)
    xmin_p,xmax_p=np.digitize([pmin[0],pmax[0]],bins)
    zmin_l,zmax_l=np.digitize([lmin[2],lmax[2]],bins)
    ymin_l,ymax_l=np.digitize([lmin[1],lmax[1]],bins)
    xmin_l,xmax_l=np.digitize([lmin[0],lmax[0]],bins)
    zmin_b,zmax_b=np.digitize([bmin[2],bmax[2]],bins)
    ymin_b,ymax_b=np.digitize([bmin[1],bmax[1]],bins)
    xmin_b,xmax_b=np.digitize([bmin[0],bmax[0]],bins)

    dc=copy.deepcopy(density)
    dc.sites=[[],[],[],[],[],[]]

    for i in xrange(zmin_p):
        for j in xrange(ymax_b):
            for k in xrange(xmax_b):
                dc.sites[1].append((k,j,i))
    for i in xrange(zmin_l):
        for j in xrange(ymin_b,ymin_p+1 ):
            for k in xrange(ymin_b,xmin_p+1):
                dc.sites[1].append((k,j,i))
    for i in xrange(zmin_l):
        for j in xrange(ymax_p,ymax_b+1 ):
            for k in xrange(ymax_p,xmax_b+1 ):
                dc.sites[1].append((k,j,i))
    for i in xrange(zmax_p,zmax_b+1):
        for j in xrange(ymax_b):
            for k in xrange(xmax_b):
                dc.sites[2].append((k,j,i))
    for i in xrange(zmax_l,zmax_b+1):
        for j in xrange(ymin_b,ymin_p+1):
            for k in xrange(ymin_b,xmin_p+1):
                dc.sites[2].append((k,j,i))
    for i in xrange(zmax_l,zmax_b+1):
        for j in xrange(ymax_p,ymax_b+1):
            for k in xrange(ymax_p,xmax_b+1):
                dc.sites[2].append((k,j,i))
    for i in xrange(zmax_l,zmax_p+1):
        for j in xrange(ymin_p,ymax_p+1):
            for k in xrange(xmin_p,xmax_p+1):
                dc.sites[3].append((k,j,i))
    for i in xrange(zmin_p,zmin_l+1):
        for j in xrange(ymin_p,ymax_p+1):
            for k in xrange(xmin_p,xmax_p+1):
                dc.sites[4].append((k,j,i))
    for i in xrange(zmin_l,zmax_l+1):
        for j in xrange(ymin_l,ymax_l+1):
            for k in xrange(xmin_l,xmax_l+1):
                dc.sites[5].append((k,j,i))
    dc.map[:,:,0:zmin_p]=1
    dc.map[0:xmin_p,0:ymin_p,0:zmin_l]=1
    dc.map[xmax_p:xmax_b,ymax_p:ymax_b,0:zmin_l]=1
    dc.map[:,:,zmax_p:zmax_b]=2
    dc.map[0:xmin_p,0:ymin_p,zmax_b:zmin_l]=2
    dc.map[xmax_p:xmax_b,ymax_p:ymax_b,zmax_b:zmin_l]=2
    dc.map[xmin_p:xmax_p,ymin_p:ymax_p,zmax_l:zmax_p]=3
    dc.map[xmin_p:xmax_p,ymin_p:ymax_p,zmin_p:zmin_l]=4
    dc.map[:,:,zmin_l:zmax_l]=5
    dc._update()
    dc._draw_map_from_sites()
    dc._annotate_sites()
    dc.export_map(labels='all')
    return dc 
