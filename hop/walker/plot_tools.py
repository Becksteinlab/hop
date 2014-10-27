import walker as wk 
import matplotlib.pyplot as plt

def fpt_distribution_plot(hopgraph,nodes,site,N):

    paths=wk.init_walker(hopgraph.graph,nodes,site,N)
    tau_list=[]
    for i in paths:
        tau=0
        for j in i:
            try:
                tau+=1/hopgraph.graph[j][i.index(j)+1]['k']    
            except:
               pass
        tau_list.append(tau)
    
    plt.plot(tau_list)
    plt.show() 
