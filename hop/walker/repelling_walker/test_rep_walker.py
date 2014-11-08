import networkx
import random
import numpy as np

def probability_updater(graph_k):
    for i in graph_k.edge:
        total_k = 0
        for j in graph_k[i]: 
            total_k += float(graph_k.edge[i][j]['k'])
        for j in graph_k[i]:
            graph_k.edge[i][j]['p']=float(graph_k.edge[i][j]['k']/total_k)
    return graph_k

def init_walker(graph,nodes,init,N):
    probability_updater(graph)
    paths=[]
    for i in xrange(N):
        csite=init
        visited_sites=[init]
        while csite not in nodes:    
            neighbors=iter(graph[csite])
            for j in neighbors:
                rand=random.random() 
                try:
                    if rand > 0 and rand < graph[csite][j]['p']:
                        csite=j
                        visited_sites.append(j)
                        break
                    else:
                        neighbors2=iter(graph[csite])
                        for i in neighbors2:
                            try:
                                if rand > graph[csite][i]['p'] and rand < graph[csite][i]['p']+graph[csite][neighbors2.next()]['p']:
                                    csite=i
                                    visited_sites.append(i)
                                    break
                                else:
                                    continue
                            except:
                                pass 
                except:
                    pass
        paths.append(visited_sites)
    return paths

def analysis(graph,nodes,init,N):
    out=init_walker(graph,nodes,init,N)
    taus=[]
    for i in out:
        tau=0
        try:
            for j in i:
                tau+=float(1/graph[j][i[i.index(j)+1]]['k'])
       #         print(tau)
        except:
            taus.append(tau)
            #print(len(taus))
            continue
    tau_avg=sum(taus)/len(taus)
    return tau_avg 
