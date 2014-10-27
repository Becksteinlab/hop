import walker as wk
import networkx as nx

three_path = nx.DiGraph()
three_path.add_edge(1,2,k=float(1))
three_path.add_edge(2,5,k=float(1))
three_path.add_edge(1,3,k=float(0.5))
three_path.add_edge(3,5,k=float(1))
three_path.add_edge(1,4,k=float(1)/float(3))
three_path.add_edge(4,5,k=float(1))

def connectivity_test():
    nodes=[5]
    init=1
    N=10
    paths=wk.init_walker(three_path,nodes,init,N)
    paths_iter=iter(paths)
    success=True
    for i in paths_iter:
        path_iter=iter(i)
        for j in i:
            try:
                if not graph.has_edge(j,j.next()):
                    success=False
                    assert success==False, "non-connected nodes"
                else:
                    continue
            except:
                continue
    
