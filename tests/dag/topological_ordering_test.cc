#include "gmock/gmock.h"

#include "kahypar/io/hypergraph_io.h"
#include "kahypar/dag/topological_ordering.h"

using ::testing::Eq;
using ::testing::ContainerEq;

namespace kahypar {
namespace dag {
static Hypergraph _loadHypergraph(const std::string& filename) {
  HypernodeID num_hypernodes;
  HyperedgeID num_hyperedges;
  HyperedgeIndexVector index_vector;
  HyperedgeVector edge_vector;
  HyperedgeWeightVector hyperedge_weights;
  HypernodeWeightVector hypernode_weights;

  io::readHypergraphFile(filename, num_hypernodes, num_hyperedges, index_vector, edge_vector,
                         &hyperedge_weights, &hypernode_weights);
  return Hypergraph(num_hypernodes, num_hyperedges, index_vector, edge_vector, 2,
                    &hyperedge_weights, &hypernode_weights);
}

static void _assertTopologicalOrderingIsTopological(
  const Hypergraph& hg,
  const std::vector<HypernodeID>& topological_ordering
) {
  ASSERT_EQ(topological_ordering.size(), hg.currentNumNodes());

  std::vector<HypernodeID> position(hg.currentNumNodes());
  for (const HypernodeID &hn : hg.nodes()) {
    position[topological_ordering[hn]] = hn;
  }

  for (const HyperedgeID &he : hg.edges()) {
    for (const HypernodeID &hh : hg.heads(he)) {
      for (const HypernodeID &ht : hg.tails(he)) {
        ASSERT_LE(position[ht], position[hh]);
      }
    }
  }
}

TEST(DAG, CyclicGraphCannotBeTopologicalOrdered) {
  Hypergraph hg = _loadHypergraph("test_instances/cyclic.hgr");
  ASSERT_THROW(calculateTopologicalOrdering(hg), CyclicGraphException);
  ASSERT_FALSE(isAcyclic(hg));
}

TEST(DAG, AcyclicGraphCanBeTopologicalOrdered) {
  Hypergraph hg = _loadHypergraph("test_instances/acyclic.hgr");
  ASSERT_NO_THROW(calculateTopologicalOrdering(hg));
  ASSERT_TRUE(isAcyclic(hg));
}

TEST(DAG, AcyclicGraphGetsTopologicalOrdered) {
  Hypergraph hg = _loadHypergraph("test_instances/acyclic.hgr");
  _assertTopologicalOrderingIsTopological(hg, calculateTopologicalOrdering(hg));
}
} // namespace dag
} // namespace kahypar