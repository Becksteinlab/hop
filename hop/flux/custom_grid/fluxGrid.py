import hop.sitemap
import MDAnalysis
import numpy as np
import copy
from custom_grid import make_grid as mkgd

def flux_graph(topology, density,N,delta=1):
    
    u=MDAnalysis.Universe(topology)
    binsx=[l for l in xrange(u.dimensions[0])]
    binsy=[l for l in xrange(u.dimensions[1])]
    binsz=[l for l in xrange(u.dimensions[2])]
    p=u.selectAtoms("protein")
    l=u.selectAtoms("resname POP")
    pmin,pmax=p.atoms.bbox()
    lmin,lmax=l.atoms.bbox()
    bmin,bmax=[0,0,0],[u.dimensions[0],u.dimensions[1],u.dimensions[2]]
    zmin_p,zmax_p=np.digitize([pmin[2],pmax[2]],binsz)
    ymin_p,ymax_p=np.digitize([pmin[1],pmax[1]],binsy)
    xmin_p,xmax_p=np.digitize([pmin[0],pmax[0]],binsx)
    zmin_l,zmax_l=np.digitize([lmin[2],lmax[2]],binsz)
    ymin_l,ymax_l=np.digitize([lmin[1],lmax[1]],binsy)
    xmin_l,xmax_l=np.digitize([lmin[0],lmax[0]],binsx)
    zmin_b,zmax_b=np.digitize([bmin[2],bmax[2]],binsz)
    ymin_b,ymax_b=np.digitize([bmin[1],bmax[1]],binsy)
    xmin_b,xmax_b=np.digitize([bmin[0],bmax[0]],binsx)

    m=mkgd(u,delta)
    grid=m[0]
    edges=m[1]
    grid_bulk=grid.copy()
    dc=copy.deepcopy(density)
    for i in xrange(N):
        """Appends N new sites to the density"""
        dc.sites.append([])

    for i in xrange(xmin_b,xmax_b):
        for j in xrange(ymin_b,ymax_b):
            for k in xrange(zmin_b,zmax_b):
                if ((i-((xmax_p+xmin_p)*0.5))**2 + (j-((ymax_p+ymin_p)*0.5))**2 + (k-zmax_l+15)**2)**(0.5) <15:
                    if k>zmax_l-15:
                        grid[i,j,k]=1
    
    for i in xrange(xmin_b,xmax_b):
        for j in xrange(ymin_b,ymax_b):
            for k in xrange(zmin_b,zmax_b):
                if ((i-((xmax_p+xmin_p)*0.5))**2 + (j-((ymax_p+ymin_p)*0.5))**2 + (k-zmin_l-11)**2)**(0.5) <15:
                    if k<zmin_l+11:
                        grid[i,j,k]=2
    
    for i in xrange(xmin_b,xmax_b):
        for j in xrange(ymin_b,ymax_b):
            for k in xrange(zmin_b,zmax_b):
                if not ((i-((xmax_p+xmin_p)*0.5))**2 + (j-((ymax_p+ymin_p)*0.5))**2 + (k-zmin_l-13)**2)**(0.5) <15:
                    if k<zmin_l:
                        grid[i,j,k]=3

    for i in xrange(xmin_b,xmax_b):
        for j in xrange(ymin_b,ymax_b):
            for k in xrange(zmin_b,zmax_b):
                if not ((i-((xmax_p+xmin_p)*0.5))**2 + (j-((ymax_p+ymin_p)*0.5))**2 + (k-zmax_l+15)**2)**(0.5) <15:
                    if k>zmax_l:
                        grid[i,j,k]=4

    for i in xrange(xmin_b,xmax_b):
        for j in xrange(ymin_b,ymax_b):
            for k in xrange(zmin_b,zmax_b):
                if k<zmax_l-16 and k>zmin_l+10:
                    grid_bulk[i,j,k]=5
    

    return [grid,edges,grid_bulk] 
