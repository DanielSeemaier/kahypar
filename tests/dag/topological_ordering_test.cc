#include "gmock/gmock.h"

#include "kahypar/dag/topological_ordering.h"

#include "dag.h"

using ::testing::Test;

namespace kahypar {
namespace dag {
class TopologicalOrderingTest : public Test {
 protected:
  void SetUp() override {
    cyclic = loadHypergraph("test_instances/cyclic.hgr");
    acyclic = loadHypergraph("test_instances/acyclic.hgr");
  }

  Hypergraph cyclic;
  Hypergraph acyclic;
};

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

TEST_F(TopologicalOrderingTest, CyclicGraphCannotBeTopologicalOrdered) {
  ASSERT_THROW(calculateTopologicalOrdering(cyclic), CyclicGraphException);
  ASSERT_FALSE(isAcyclic(cyclic));
}

TEST_F(TopologicalOrderingTest, AcyclicGraphCanBeTopologicalOrdered) {
  ASSERT_NO_THROW(calculateTopologicalOrdering(acyclic));
  ASSERT_TRUE(isAcyclic(acyclic));
}

TEST_F(TopologicalOrderingTest, AcyclicGraphGetsTopologicalOrdered) {
  _assertTopologicalOrderingIsTopological(acyclic, calculateTopologicalOrdering(acyclic));
}
} // namespace dag
} // namespace kahypar