import hop.graph
import numpy as np
import utilities

def get_matrix(hopgraph,quantity='p'):
    utilities.graph_updater(hopgraph.graph)
    transition_matrix=np.zeros((len(hopgraph.graph.nodes()),len(hopgraph.graph.nodes())),dtype=np.float64)
    for i in hopgraph.graph.nodes():
        for j in hopgraph.graph.nodes():
            if hopgraph.graph.has_edge(i,j):
                i_index=hopgraph.graph.nodes().index(i)
                j_index=hopgraph.graph.nodes().index(j)
                transition_matrix.itemset((i_index,j_index), hopgraph.graph[i][j][quantity])
    return transition_matrix

def mfpt_eq(transition_matrix, tau):

    T=transition_matrix
    T_tau=T*tau

    a=np.zeros(T.shape) 
    b=np.zeros(T.shape)
    T_sum_matrix=np.zeros(T.shape,dtype=np.float64)
    for i in xrange(T.shape[0]):
        for j in xrange(T.shape[0]):
            T_sum=0
            for k in xrange(T.shape[0]):
                if k!=j:
                    T_sum+=-T.item((i,k))*tau
            T_sum_matrix.itemset((i,j),T_sum)
    for i in xrange(T.shape[0]):
        for j in xrange(T.shape[0]):
            a.itemset((i,j),1-T_sum_matrix.item((i,j)))
            b.itemset((i,j), T.item((i,j))*tau-T_sum_matrix.item((i,j)))     
    s=np.zeros(T.shape)
    for i in xrange(T.shape[0]):
        s[i]=np.linalg.solve(a,b[i])
    return [s,T_sum_matrix,a,b]
