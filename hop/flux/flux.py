import hop.sitemap
import MDAnalysis
import numpy as np
import copy

def flux_graph(topology, density,N):
    
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

    dc=copy.deepcopy(density)
 #   del dc.sites
#    del dc.map
#    dc.sites=[[],[],[],[],[],[]]
    for i in xrange(N):
        """Appends N new sites to the density"""
        dc.sites.append([])

    for i in xrange(xmin_b,xmax_b):
        for j in xrange(ymin_b,ymax_b):
            for k in xrange(zmin_b,zmax_b):
                if ((i-((xmax_p+xmin_p)*0.5))**2 + (j-((ymax_p+ymin_p)*0.5))**2 + (k-zmax_l+15)**2)**(0.5) <15:
                    dc.sites[-4].append((i,j,k))
    
    for i in xrange(xmin_b,xmax_b):
        for j in xrange(ymin_b,ymax_b):
            for k in xrange(zmin_b,zmax_b):
                if ((i-((xmax_p+xmin_p)*0.5))**2 + (j-((ymax_p+ymin_p)*0.5))**2 + (k-zmin_l-11)**2)**(0.5) <15:
                    dc.sites[-3].append((i,j,k))
    
    for i in xrange(xmin_b,xmax_b):
        for j in xrange(ymin_b,ymax_b):
            for k in xrange(zmin_b,zmax_b):
                if not ((i-((xmax_p+xmin_p)*0.5))**2 + (j-((ymax_p+ymin_p)*0.5))**2 + (k-zmin_l-13)**2)**(0.5) <15:
                    if k<zmin_l:
                        dc.sites[-2].append((i,j,k))

    for i in xrange(xmin_b,xmax_b):
        for j in xrange(ymin_b,ymax_b):
            for k in xrange(zmin_b,zmax_b):
                if not ((i-((xmax_p+xmin_p)*0.5))**2 + (j-((ymax_p+ymin_p)*0.5))**2 + (k-zmax_l+15)**2)**(0.5) <15:
                    if k>zmax_l:
                        dc.sites[-1].append((i,j,k))

#    for i in xrange(xmin_b,xmax_b):
 #       for j in xrange(ymin_b,ymax_b):
  #          for k in xrange(zmin_b,zmax_b):
   #             c=(i,j,k)
    #            if c not in set(dc.sites[1]) or c not in set(dc.sites[2]) or c not in set(dc.sites[3]) or c not in set(dc.sites[4]):
        #                dc.sites[5].append(c)
    """
    for i in xrange(xmin_b,xmax_b):
        for j in xrange(ymin_b,ymax_b):
            for k in xrange(zmin_b,zmax_b):
                if not ((i-((xmax_p+xmin_p)*0.5))**2 + (j-((ymax_p+ymin_p)*0.5))**2 + (k-zmin_l-13)**2)**(0.5) <15:
                    if ((i-(xmax_b+xmin_b*0.5))**2 + (j-(ymax_b+ymin_b*0.5))**2+(k-(zmax_b+zmin_b*0.5))**2)**(0.5) < zmax_b/2 and k>zmax_b/2  :
                    dc.sites[3].append((i,j,k))
                    """
    """
    dc.map[0:xmin_p,0:ymin_p,0:zmin_l]=1
    dc.map[xmax_p:xmax_b,ymax_p:ymax_b,0:zmin_l]=1
    dc.map[:,:,zmax_p:zmax_b]=2
    dc.map[0:xmin_p,0:ymin_p,zmax_b:zmin_l]=2
    dc.map[xmax_p:xmax_b,ymax_p:ymax_b,zmax_b:zmin_l]=2
    dc.map[xmin_p:xmax_p,ymin_p:ymax_p,zmax_l:zmax_p]=3
    dc.map[xmin_p:xmax_p,ymin_p:ymax_p,zmin_p:zmin_l]=4
    dc.map[:,:,zmin_l:zmax_l]=5"
         """
    for i in dc.sites[-4]:
        dc.map[i[0],i[1],i[2]]=1000
    for i in dc.sites[-3]:
        dc.map[i[0],i[1],i[2]]=2
    dc._update()
    dc._draw_map_from_sites()
    dc._annotate_sites()
    dc.export_map(labels=[281])
    return dc 
