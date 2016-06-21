# $Id$
# stuff used during debugging/testing


def load_my_sims():
    """Dict of CombinedGraph instances as input for HeatmapAnalysis.
    Hard coded loader for testing.
    """
    ligand = {0:'apo',1:'plm'}  # state --> ligand
    hopgraphs = {}

    # hard coded list for testing, use expressive keys
    # Single hop graphs
    hopgraphs['PME_s_apo'] = hop.graph.HoppingGraph(filename='/Users/oliver/Biop/Projects/WaterNetworks/testcases/1IFC/hopgraph.pickle')
    hopgraphs['PME_s_plm'] = hop.graph.HoppingGraph(filename='/Users/oliver/Biop/Projects/WaterNetworks/testcases/2IFB/hopgraph.pickle')

    # combined hopgraphs. G0 == apo, G1 == holo==plm
    cg = {}
    cg['R13_1'] = hop.graph.CombinedGraph(filename='GSBP/R13/analysis/cg_1IFC_2IFB_1.pickle')
    cg['R15_0'] = hop.graph.CombinedGraph(filename='GSBP/R15/analysis/cg_1IFC_2IFB_0.pickle')
    cg['R15_1'] = hop.graph.CombinedGraph(filename='GSBP/R15/analysis/cg_1IFC_2IFB_1.pickle')
    cg['PME_0'] = hop.graph.CombinedGraph(filename='LAMMPS/100mM/analysis/cg_1IFC_2IFB.pickle')
    cg['PME_TAP']=hop.graph.CombinedGraph(filename='LAMMPS/100mM/analysis/cg_TAP_1IFC_2IFB.pickle')
    cg['PME_CRBP'] = hop.graph.CombinedGraph(filename='../../../CRBPII/comparison/analysis/cg_1OPA_1OPB.pickle')
    
    for sim,combgraph in cg.items():
        for state,hopgraph in enumerate(combgraph.graphs):  # assuming a CombinedGraph
            #print "sim = %(sim)s  state = %(state)d"  % locals()
            sim_id = str(sim)+'_'+ligand[state]
            hopgraphs[sim_id] = hopgraph

    return hopgraphs

