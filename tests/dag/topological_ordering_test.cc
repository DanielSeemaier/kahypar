#include "gmock/gmock.h"

#include "kahypar/dag/topological_ordering.h"

#include "dag.h"

using ::testing::Eq;
using ::testing::ContainerEq;

namespace kahypar {
namespace dag {
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