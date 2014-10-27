import hop.graph

def symmetry_test(hopgraph=None,graph=None):
    diff_list=[]
    non_symmetric=[]
    if hopgraph !=None:
        for i in hopgraph.graph.edges():
            if i[0] in hopgraph.graph[i[1]]: 
                diff=hopgraph.graph[i[0]][i[1]]['k']-hopgraph.graph[i[1]][i[0]]['k']
                diff_list.append(diff)
            else:
                non_symmetric.append(i)
    avg=sum(diff_list)/len(diff_list)
    return [avg,non_symmetric]
