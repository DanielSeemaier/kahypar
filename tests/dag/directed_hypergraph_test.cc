#include "gmock/gmock.h"

#include "kahypar/definitions.h"
#include "kahypar/partition/context.h"
#include "kahypar/dag/topological_ordering.h"
#include "kahypar/io/partitioning_output.h"

#include "dag.h"

using ::testing::Eq;

namespace kahypar {
namespace dag{
TEST(DirectedHypergraphTest, DirectedHypergraphIsDirected) {
  // c17 is a directed hypergraph where each edge has one head and two tails
  auto hgr = loadHypergraph("test_instances/c17.hgr");
  ASSERT_THAT(hgr.isDirected(), Eq(true));
  for (const HyperedgeID& he : hgr.edges()) {
    ASSERT_THAT(hgr.edgeNumHeads(he), Eq(1));
    ASSERT_THAT(hgr.edgeNumTails(he), Eq(2));
  }
  for (const auto& hn : {0, 1, 2, 3, 5, 9}) {
    ASSERT_THAT(hgr.nodeNumHeads(hn), Eq(1));
  }
  for (const auto& hn : {1, 2, 10}) {
    ASSERT_THAT(hgr.nodeNumTails(hn), Eq(2));
  }
  for (const auto& hn : {0, 4, 5, 6, 7, 8}) {
    ASSERT_THAT(hgr.nodeNumTails(hn), Eq(1));
  }
}

TEST(DirectedHypergraphTest, UndirectedHypergraphIsUndirected) {
  auto hgr = loadHypergraph("test_instances/undirected.hgr");
  ASSERT_THAT(hgr.isDirected(), Eq(false));
}
} // namespace dag
} // namespace kahypar

