import hop.graph
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

def rate_distribution(hopgraph,bins=10,alpha=0.5,log=False):
    
    h=hopgraph
    h.filter(exclude={'bulk':True,'unconnected':True,'outliers':True})
    rates=[]
    for i in h.filtered_graph:
        for j in h.filtered_graph[i]:
            rates.append(h.filtered_graph[i][j]['k'])
    
    plt.hist(rates,bins,alpha=alpha,log=log)
