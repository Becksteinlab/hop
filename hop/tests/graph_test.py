#!/usr/bin/env python
"""%prog [options] -g ORIGINAL_GRAPH -t COMPARISON_GRAPH

Compare the attributes of two hopgraph objects to see if they are identical. Used for integration purposes.


"""

import os.path, errno
import numpy 
import MDAnalysis
import hop.interactive
import hop.graph
from hop.utilities import unlink_f, mkdir_p

import logging
logger = logging.getLogger('MDAnalysis.app')


if __name__ == "__main__":
    import sys
    from optparse import OptionParser

    parser = OptionParser(usage=__doc__)
    parser.add_option("-g", "--graph", dest="graph",
                      metavar="FILE",
                      help="Graph file (e.g. hopgraph.pickle) for analysis [%default]")
    parser.add_option("-D", "--analysisdir", dest="analysisdir",
                      metavar="DIRNAME",
                      help="results will be stored under DIRNAME/(basename DIR)  [%default]")
    parser.add_option("-t", "--analysisgraph", dest="analysisgraph",
                      metavar="COMPFILE",
                      help="Comparison graph file  [%default]")

    opts,args = parser.parse_args()

    analysisdir = os.path.join(opts.analysisdir)
    
    new_graph = hop.graph.HoppingGraph(filename=opts.graph)
    old_graph = hop.graph.HoppingGraph(filename=opts.analysisgraph)

    new_graph.filter(exclude={'unconnected':True,'bulk':True,'outliers':True})
    old_graph.filter(exclude={'unconnected':True,'bulk':True,'outliers':True})
    new_graph.compute_site_occupancy()
    new_graph.compute_site_times
    old_graph.compute_site_occupancy()
    old_graph.compute_site_times

    MDAnalysis.start_logging()

    if len(args) != 0:
        logger.fatal("This command only accepts option arguments. See --help.")
        sys.exit(1)

    graph = os.path.abspath(opts.graph)
    if not os.path.exists(graph):
        errmsg = "Graph %(Graph)r not found; (use --graph)" % vars()
        logger.fatal(errmsg)
        raise IOError(errno.ENOENT, errmsg)

    analysisgraph = os.path.abspath(opts.analysisgraph)
    if not os.path.exists(analysisgraph):
        errmsg = "Analysis graph %(analysisgraph)r not found; (use -t)" % vars()
        logger.fatal(errmsg)
        raise IOError(errno.ENOENT, errmsg)

    logger.debug("graph   = %(graph)r", vars())

    startdirectory = os.path.abspath(os.path.curdir)
     
    numpy.testing.assert_array_almost_equal(new_graph.occupancy_avg,old_graph.occupancy_avg,err_msg="Occupancies are unequal")

    numpy.testing.assert_array_almost_equal(new_graph.filtered_graph.edges(),old_graph.filtered_graph.edges(),err_msg="Edges are unequal")

    numpy.testing.assert_array_almost_equal(new_graph.filtered_graph.nodes(),old_graph.filtered_graph.nodes(),err_msg="Node lists are unequal")

    if new_graph.filtered_graph.degree() != old_graph.filtered_graph.degree():
        print('Unequal Degrees')

    if new_graph.filtered_graph.number_of_selfloops() != old_graph.filtered_graph.number_of_selfloops():
        print('Unequal number of self loops')

    if new_graph.filtered_graph.order() != old_graph.filtered_graph.order():
        print('Unequal graph order')
    
    if new_graph.filtered_graph.size() != old_graph.filtered_graph.size():
        print('Unequal graph sizes')

MDAnalysis.stop_logging()

