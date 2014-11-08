import hop.sitemap
import copy
import MDAnalysis as MDA

def membrane_graph_gen(sample_density,universe):
    
    dc=copy.deepcopy(sample_density)
    for i in sample_density.sites[1]:
        try:
            del dc.sites[1][1]
        except:
            pass

    u=MDA.Universe(universe)
    bmin, bmax = u.atoms.bbox()
    box=bmax - bmin

    dim=u.dimensions
    bottom=dim[2]/2 - box[2]/2

    for i in xrange(dim[0]):
        for j in xrange(dim[1]):
            for k in [l for l in xrange(dim[2]) if l>box[2]]:
                dc.sites[1].append((i,j,k))
            for k in [l for l in xrange(dim[2]) if l<bottom]:
                dc.sites[1].append((i,j,k))    
    return dc    
