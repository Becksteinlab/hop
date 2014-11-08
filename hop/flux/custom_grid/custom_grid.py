import hop.sitemap
import hop.density
from hop.sitemap import fixedwidth_bins
import MDAnalysis
import numpy as np

def make_grid(universe,delta):

    u=universe    
    BINS=fixedwidth_bins(delta,np.min(u.atoms.coordinates(),axis=0),np.max(u.atoms.coordinates(),axis=0))

    arange=zip(BINS['min'],BINS['max'])

    bins=BINS['Nbins']

    grid,edges = np.histogramdd(np.zeros((1,3)),bins=bins,range=arange,normed=False)
    grid *= 0.0
    return [grid,edges]
