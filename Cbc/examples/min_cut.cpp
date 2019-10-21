/**
  The original version of this file is taken from the MiniZinc source code
  www.minizinc.org
  */

/*
 *  Main authors:
 *     Gleb Belov <gleb.belov@monash.edu>
 */

#include <stdexcept>

#include "min_cut.h"

#define COMPILE_BOOST_MINCUT
#ifdef COMPILE_BOOST_MINCUT
  #include <boost/graph/adjacency_list.hpp>
  #include <boost/graph/graph_traits.hpp>
  #include <boost/graph/one_bit_color_map.hpp>
  #include <boost/graph/stoer_wagner_min_cut.hpp>
  #include <boost/property_map/property_map.hpp>
  #include <boost/typeof/typeof.hpp>
#endif

void Algorithms::MinCut::solve() {
#ifndef COMPILE_BOOST_MINCUT
  throw std::runtime_error( "MIP/circuit: Subtour Elimination Constraints: MinCut::solve not compiled. "
  "Compile with COMPILE_BOOST_MINCUT (needs boost), or use -D nSECcuts=0 for flattening." );
#else
  typedef std::pair< unsigned long, unsigned long > edge_t;

  // A graphic of the min-cut is available at <http://www.boost.org/doc/libs/release/libs/graph/doc/stoer_wagner_imgs/stoer_wagner.cpp.gif>
  using namespace std;
  
  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
    boost::no_property, boost::property<boost::edge_weight_t, double> > undirected_graph;
  typedef boost::property_map<undirected_graph, boost::edge_weight_t>::type weight_map_type;
  typedef boost::property_traits<weight_map_type>::value_type weight_type;
  
  // define the 16 edges of the graph. {3, 4} means an undirected edge between vertices 3 and 4.
//   edge_t edges[] = {{3, 4}, {3, 6}, {3, 5}, {0, 4}, {0, 1}, {0, 6}, {0, 7},
//     {0, 5}, {0, 2}, {4, 1}, {1, 6}, {1, 5}, {6, 7}, {7, 5}, {5, 2}, {3, 4}};
  std::vector< edge_t > edges_array;  // = edges;
//   std::vector< weight_type > ws_array;
  
  edges_array.reserve( edges.size() );
  for ( const auto& edg: edges )
    edges_array.push_back( std::pair<unsigned long, unsigned long>( edg.first, edg.second ) );
  
  // for each of the 16 edges, define the associated edge weight. ws[i] is the weight for the edge
  // that is described by edges[i].
//   weight_type ws[] = {0.2, 3.1, 1.3, 3.7, 1.5, 2.6, 6.44, 1.26, 8.77, 1.29, 1.95, 80.74, 2.23, 1.94, 1.23, 4.957};
  
  // construct the graph object. 8 is the number of vertices, which are numbered from 0
  // through 7, and 16 is the number of edges.
  undirected_graph g(edges_array.data(), edges_array.data() + edges_array.size(),
                     weights.data(), nNodes, edges_array.size());
  
  // define a property map, `parities`, that will store a boolean value for each vertex.
  // Vertices that have the same parity after `stoer_wagner_min_cut` runs are on the same side of the min-cut.
  BOOST_AUTO(parities00, boost::make_one_bit_color_map(num_vertices(g), get(boost::vertex_index, g)));
  
  // run the Stoer-Wagner algorithm to obtain the min-cut weight. `parities` is also filled in.
  wMinCut = boost::stoer_wagner_min_cut(g, get(boost::edge_weight, g), boost::parity_map(parities00));
  
  assert( nNodes==num_vertices(g) );
  
//   parities = parities00;
  parities.clear();
  parities.reserve( nNodes );
  for (size_t i = 0; i < num_vertices(g); ++i) {
    parities.push_back (get(parities00, i));
  }
#endif
}
