#include <algorithm>

#include "gmock/gmock.h"

#include "dag.h"

#include "kahypar/partition/refinement/policies/fm_stop_policy.h"
#include "kahypar/dag/fix_acyclicity.h"

using ::testing::TestWithParam;
using ::testing::Values;
using ::testing::Ge;
using ::testing::Eq;

namespace kahypar {
namespace dag {
TEST(FixAcyclicityTest, DetectsCycleOfLengthTwo) {
  auto hg = loadHypergraph("test_instances/c17.hgr");
  hg.changeK(2);
  hg.setNodePart(0, 0);
  hg.setNodePart(1, 1);
  hg.setNodePart(2, 1);
  hg.setNodePart(3, 1);
  hg.setNodePart(4, 1);
  hg.setNodePart(5, 1);
  hg.setNodePart(6, 0);
  hg.setNodePart(7, 0);
  hg.setNodePart(8, 0);
  hg.setNodePart(9, 0);
  hg.setNodePart(10, 1);

  Context context;
  context.partition.k = 2;
  CyclicQuotientGraph cqg(hg, context);

  ASSERT_FALSE(AdjacencyMatrixQuotientGraph<DFSCycleDetector>(hg, context).isAcyclic());
  const auto cycle = internal::findCycle(cqg);
  ASSERT_THAT(cycle.size(), Eq(2));
  ASSERT_THAT(cycle.front().second, Eq(cycle.back().first));
  ASSERT_THAT(cycle.front().first, Eq(cycle.back().second));
}

TEST(FixAcyclicityTest, DetectsCycleOfLengthThree) {
  auto hg = loadHypergraph("test_instances/c17.hgr");
  hg.changeK(4);
  hg.setNodePart(0, 1);
  hg.setNodePart(1, 3);
  hg.setNodePart(2, 3);
  hg.setNodePart(3, 3);
  hg.setNodePart(4, 2);
  hg.setNodePart(5, 3);
  hg.setNodePart(6, 1);
  hg.setNodePart(7, 0);
  hg.setNodePart(8, 2);
  hg.setNodePart(9, 2);
  hg.setNodePart(10, 2);

  Context context;
  context.partition.k = 4;
  CyclicQuotientGraph cqg(hg, context);

  ASSERT_FALSE(AdjacencyMatrixQuotientGraph<DFSCycleDetector>(hg, context).isAcyclic());
  const auto cycle = internal::findCycle(cqg);
  ASSERT_THAT(cycle.size(), Eq(3));
  ASSERT_THAT(cycle[0].second, Eq(cycle[1].first));
  ASSERT_THAT(cycle[1].second, Eq(cycle[2].first));
  ASSERT_THAT(cycle[2].second, Eq(cycle[0].first));

  const auto chosen_edge = internal::pickEdge(cqg);
  ASSERT_TRUE(std::find(cycle.begin(), cycle.end(), chosen_edge) != cycle.end());
  for (const auto& edge : cycle) {
    ASSERT_THAT(cqg.edgeWeight(edge.first, edge.second),
      Ge(cqg.edgeWeight(chosen_edge.first, chosen_edge.second)));
  }
}
} // namespace dag
} // namespace kahypar