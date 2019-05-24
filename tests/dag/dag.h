#pragma once

#include "kahypar/io/hypergraph_io.h"

namespace kahypar {
namespace dag {
static Hypergraph _loadHypergraph(const std::string& filename, const bool directed = true,
                                  const HypernodeID num_heads_per_hyperedge = 1) {
  HypernodeID num_hypernodes;
  HyperedgeID num_hyperedges;
  HyperedgeIndexVector index_vector;
  HyperedgeVector edge_vector;
  HyperedgeWeightVector hyperedge_weights;
  HypernodeWeightVector hypernode_weights;

  io::readHypergraphFile(filename, num_hypernodes, num_hyperedges, index_vector, edge_vector,
                         &hyperedge_weights, &hypernode_weights);
  return Hypergraph(num_hypernodes, num_hyperedges, index_vector, edge_vector, 2,
                    &hyperedge_weights, &hypernode_weights, directed, num_heads_per_hyperedge);
}
} // namespace dag
} // namespace kahypar