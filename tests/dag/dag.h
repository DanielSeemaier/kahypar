#pragma once

#include <algorithm>

#include "kahypar/datastructure/hypergraph.h"
#include "kahypar/io/hypergraph_io.h"

namespace kahypar {
namespace dag {
static Hypergraph loadHypergraph(const std::string& filename) {
  HypernodeID num_hypernodes;
  HyperedgeID num_hyperedges;
  HyperedgeIndexVector index_vector;
  HyperedgeVector edge_vector;
  HyperedgeWeightVector hyperedge_weights;
  HypernodeWeightVector hypernode_weights;
  NumHeadsVector num_heads_vector;
  bool is_directed;

  io::readHypergraphFile(filename, num_hypernodes, num_hyperedges, index_vector, edge_vector,
                         is_directed, num_heads_vector, &hyperedge_weights, &hypernode_weights);
  return Hypergraph(num_hypernodes, num_hyperedges, index_vector, edge_vector, is_directed, num_heads_vector, 2,
                    &hyperedge_weights, &hypernode_weights);
}

void placeAllHypernodesInPartition(Hypergraph& hg, const PartitionID part = 0) {
  for (const HypernodeID& hn : hg.nodes()) {
    hg.setNodePart(hn, part);
  }
}

template<typename Iterator>
auto toVec(const std::pair<Iterator, Iterator>&& iter_pair) -> std::vector<typename Iterator::value_type> {
  std::vector<typename Iterator::value_type> values;
  values.insert(values.begin(), iter_pair.first, iter_pair.second);
  return values;
}
} // namespace dag
} // namespace kahypar