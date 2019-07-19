#pragma once

#include <algorithm>

#include "kahypar/dag/topological_ordering.h"
#include "kahypar/datastructure/hypergraph.h"
#include "kahypar/io/hypergraph_io.h"
#include "kahypar/partition/context.h"
#include "kahypar/utils/randomize.h"

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

class BaseDAGTest {
 protected:
  virtual void loadGraph(const std::string& filename, const PartitionID k) {
    hg = loadHypergraph(filename);
    partitionUsingTopologicalOrdering(k);
    context.partition.k = k;
    context.partition.epsilon = 0.03;
    context.setupPartWeights(hg.totalWeight());
  }

  void partitionUsingTopologicalOrdering(const PartitionID k) {
    Randomize::instance().setSeed(0);
    auto ordering = calculateTopologicalOrdering(hg);
    PartitionID part = 0;
    HypernodeID nodes_per_part = hg.initialNumNodes() / k + 1;
    HypernodeID nodes_in_cur_part = 0;

    hg.resetPartitioning();
    hg.changeK(k);
    for (const HypernodeID& hn : ordering) {
      hg.setNodePart(hn, part);
      ++nodes_in_cur_part;
      if (nodes_in_cur_part == nodes_per_part) {
        ++part;
        nodes_in_cur_part = 0;
      }
    }
  }

  std::vector<Hypergraph::Memento> contractArbitrarily(const std::size_t iterations_limit,
                                                       const std::size_t contraction_limit,
                                                       const bool keep_acyclic = false) {
    std::vector<Hypergraph::Memento> contractions;
    auto& rand = Randomize::instance();
    std::vector<bool> matched(hg.initialNumNodes());
    std::vector<HypernodeID> current_hns;
    for (const HypernodeID& hn : hg.nodes()) {
      current_hns.push_back(hn);
    }

    std::size_t num_contraction = 0;
    for (std::size_t it = 0; it < iterations_limit; ++it) {
      for (const HypernodeID& hn : current_hns) {
        if (matched[hn]) {
          continue;
        }

        for (const HyperedgeID& he : hg.incidentEdges(hn)) {
          for (const HypernodeID& pin : hg.pins(he)) {
            if (!matched[pin] && pin != hn && hg.partID(pin) == hg.partID(hn)) {
              matched[hn] = true;
              matched[pin] = true;
              contractions.push_back(hg.contract(hn, pin));
              break;
            }
          }
          if (matched[hn]) {
            break;
          }
        }
      }

      ++num_contraction;
      if (num_contraction == contraction_limit) {
        break;
      }

      current_hns.clear();
      for (const HypernodeID& hn : hg.nodes()) {
        current_hns.push_back(hn);
      }
      matched.clear();
      matched.resize(hg.initialNumNodes(), false);
    }

    hg.initializeNumCutHyperedges();
    return contractions;
  }

  Hypergraph hg;
  Context context;
};
} // namespace dag
} // namespace kahypar