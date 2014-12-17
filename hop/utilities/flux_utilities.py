import hop.graph
import numpy as np
import random
import networkx as nx

def init_prob_vect(hgraph):
    prob_vect_dict={}
    for i in hgraph:    
        tot_prob=0
        for j in hgraph[i]:
            tot_prob+=hgraph[i][j]['p']
        for k in hgraph[i]:
            prob_vect=[]
            for s in xrange(len(hgraph[i])):
                for n in xrange(len(hgraph[i])):
                    numerator=0
                    while n <s:
                        numerator+=hgraph[i][k]['p']
                        print(numerator)
                        prob_vect.append(numerator/tot_prob)
                        print(numerator/tot_prob)
                    prob_vect_dict[i]=prob_vect
                    print(prob_vect_dict[i])
    return prob_vect_dict 

def init_graph(hopgraph,do_filter=False):

    if do_filter==True:

        hopgraph.filter(exclude={"outliers":True,"bulk":True,"unconnected":True})
        f=hopgraph.filtered_graph
        return f
    else:
        graph=graph_updater(hopgraph)
        return graph

def init_state(nodes,occupancy=1,hopgraph=None,graph=None):
    if hopgraph !=None:
        state=np.zeros(len(hopgraph.graph.nodes()))
    elif graph !=None:
        state=np.zeros(len(graph.nodes()))
    for i in nodes:
        state.itemset(i-1,occupancy)
    return state

def graph_updater(graph_k):

    for i in graph_k.edge:
        total_k = 0
        for j in graph_k[i]: 
            total_k += float(graph_k.edge[i][j]['k'])
        for j in graph_k[i]:
            graph_k.edge[i][j]['p']=float(graph_k.edge[i][j]['k']/total_k)

    for i in graph_k.edge:
        total_k = 0
        for j in graph_k[i]:
            graph_k.edge[i][j]['t']=1/float(graph_k.edge[i][j]['k'])
    return graph_k

def init_test_graph(nodenumber=None,edge_weights=1,uniform=True,random_graph=True,interval=(0,1),fully_connected=True):
    G=nx.DiGraph()
    if fully_connected is True:
        for i in xrange(nodenumber):
            G.add_node(i+1)
    if uniform is True:
        for i in G.nodes():
            G.add_edge(i,i+1,k=1)
    if random is True and uniform is True:
        for i in G.nodes():
            rand=random.uniform(interval[0],interval[1])
            G.add_edge(i,i+1,k=rand)
    if random_graph is True and uniform is False:
        for j in xrange(random.randint(1,len(G.nodes()))):
            rand=random.uniform(interval[0],interval[1])
            G.add_edge(random.choice(G.nodes()),random.choice(G.nodes()),k=rand)
    return G    
         
    
