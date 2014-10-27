# coding: utf-8
import hop.graph
import walker as wk
C= hop.graph.HoppingGraph(filename='../../../../../Projects/Mhp1/out_open/hop_1.7/hopgraph.pickle')
C.filter(exclude={'bulk':True,'outliers':True,'unconnected':True})
nodes=[k for k in C.graph.nodes() if C.site_properties[k][4][2] <= 30]