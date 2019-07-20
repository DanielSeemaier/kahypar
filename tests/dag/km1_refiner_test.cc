#include "gmock/gmock.h"

#include "kahypar/definitions.h"
#include "kahypar/partition/context.h"
#include "kahypar/dag/topological_ordering.h"
#include "kahypar/io/partitioning_output.h"
#include "kahypar/partition/refinement/acyclic_km1_refiner.h"

#include "dag.h"

using ::testing::Eq;
using ::testing::Ne;
using ::testing::Test;
using ::testing::TestWithParam;
using ::testing::UnorderedElementsAre;
using ::testing::IsEmpty;
using ::testing::Values;
using ::testing::_;
using ::testing::Le;

namespace kahypar {
namespace dag {
class KMinusOneRefinerTest : public BaseDAGTest, public TestWithParam<const char*> {
 protected:
  void SetUp() override {
    std::string filename;
    PartitionID k;
    std::tie(filename, k) = parseParams();
    loadGraph(filename, k);
  }

  void loadGraph(const std::string& filename, const PartitionID k) override {
    BaseDAGTest::loadGraph(filename, k);
    refiner = std::make_unique<AcyclicKMinusOneRefiner>(hg, context);
  }

  void runRefiner(const double target_epsilon, std::vector<HypernodeID> refinement_nodes = {}) {
    context.partition.epsilon = target_epsilon;
    context.setupPartWeights(hg.totalWeight());
    if (refinement_nodes.empty()) {
      refinement_nodes = borderNodes();
    }
    Metrics metrics{metrics::hyperedgeCut(hg), metrics::km1(hg), metrics::imbalance(hg, context)};
    refiner->refine(refinement_nodes, {0, 0}, {}, metrics);

    ASSERT_THAT(metrics::km1(hg), Eq(metrics.km1));
    ASSERT_THAT(metrics::imbalance(hg, context), Eq(metrics.imbalance));
    ASSERT_THAT(metrics::imbalance(hg, context), Le(target_epsilon));
  }

 private:
  std::pair<std::string, PartitionID> parseParams() const {
    const std::string param = GetParam();
    std::string graph_filename;
    PartitionID k;
    std::stringstream(param) >> graph_filename >> k;
    return {graph_filename, k};
  }

 protected:
  std::unique_ptr<AcyclicKMinusOneRefiner> refiner;
};

TEST_P(KMinusOneRefinerTest, CanRefineBalancedPartition) {
  refiner->initialize(0);
  runRefiner(context.partition.epsilon, borderNodes());
}

TEST_P(KMinusOneRefinerTest, CanRefineUnbalancedPartition) {
  constexpr double FINAL_EPSILON = 0.03;

  applyImbalancedPartition(std::max(context.partition.k / 4, 1), 0.75);
  refiner->initialize(0);
  runRefiner(FINAL_EPSILON);
  ASSERT_THAT(metrics::imbalance(hg, context), Le(FINAL_EPSILON));
}

TEST_P(KMinusOneRefinerTest, CanRefineBalancedPartitionDuringUncoarsening) {
  auto contractions = contractArbitrarily(3, 250, true);
  refiner->initialize(0);

  constexpr double FINAL_EPSILON = 0.03;
  context.partition.epsilon = FINAL_EPSILON;
  context.partition.final_epsilon = FINAL_EPSILON;
  context.setupPartWeights(hg.totalWeight());
  Metrics metrics{metrics::hyperedgeCut(hg), metrics::km1(hg), metrics::imbalance(hg, context)};

  for (auto rit = contractions.crbegin(); rit != contractions.crend(); ++rit) {
    refiner->preUncontraction(rit->u);
    hg.uncontract(*rit);
    refiner->postUncontraction(rit->u, {rit->v});
    std::vector<HypernodeID> refinement_nodes{rit->u, rit->v};
    refiner->refine(refinement_nodes, {0, 0}, {}, metrics);

    ASSERT_THAT(metrics::km1(hg), Eq(metrics.km1));
    ASSERT_TRUE(AdjacencyMatrixQuotientGraph<DFSCycleDetector>(hg, context).isAcyclic());
  }

  ASSERT_THAT(metrics::imbalance(hg, context), Le(FINAL_EPSILON));
}

TEST_P(KMinusOneRefinerTest, CanRefineImbalancedPartitionDuringUncoarsening) {
  applyImbalancedPartition(std::max(context.partition.k / 4, 1), 0.75);
  auto contractions = contractArbitrarily(3, 250, true);
  refiner->initialize(0);

  constexpr double FINAL_EPSILON = 0.03;
  context.partition.epsilon = 1.0;
  context.partition.final_epsilon = FINAL_EPSILON;
  context.setupPartWeights(hg.totalWeight());
  Metrics metrics{metrics::hyperedgeCut(hg), metrics::km1(hg), metrics::imbalance(hg, context)};

  for (auto rit = contractions.crbegin(); rit != contractions.crend(); ++rit) {
    refiner->preUncontraction(rit->u);
    hg.uncontract(*rit);
    refiner->postUncontraction(rit->u, {rit->v});
    std::vector<HypernodeID> refinement_nodes{rit->u, rit->v};
    refiner->refine(refinement_nodes, {0, 0}, {}, metrics);

    ASSERT_THAT(metrics::km1(hg), Eq(metrics.km1));
    ASSERT_TRUE(AdjacencyMatrixQuotientGraph<DFSCycleDetector>(hg, context).isAcyclic());
  }

  refiner->printSummary();
  ASSERT_THAT(metrics::imbalance(hg, context), Le(FINAL_EPSILON));
}

INSTANTIATE_TEST_CASE_P(GRAPH_C17_K_2, KMinusOneRefinerTest, Values("test_instances/c17.hgr 2"));

INSTANTIATE_TEST_CASE_P(GRAPH_C17_K_4, KMinusOneRefinerTest, Values("test_instances/c17.hgr 4"));

INSTANTIATE_TEST_CASE_P(GRAPH_C3540_K_2, KMinusOneRefinerTest, Values("test_instances/c3540.hgr 2"));

INSTANTIATE_TEST_CASE_P(GRAPH_C3540_K_4, KMinusOneRefinerTest, Values("test_instances/c3540.hgr 4"));

INSTANTIATE_TEST_CASE_P(GRAPH_C3540_K_64, KMinusOneRefinerTest, Values("test_instances/c3540.hgr 64"));

INSTANTIATE_TEST_CASE_P(GRAPH_C7552_K_2, KMinusOneRefinerTest, Values("test_instances/c7552.hgr 2"));

INSTANTIATE_TEST_CASE_P(GRAPH_C7552_K_4, KMinusOneRefinerTest, Values("test_instances/c7552.hgr 4"));

INSTANTIATE_TEST_CASE_P(GRAPH_C7552_K_8, KMinusOneRefinerTest, Values("test_instances/c7552.hgr 8"));

INSTANTIATE_TEST_CASE_P(GRAPH_C7552_K_64, KMinusOneRefinerTest, Values("test_instances/c7552.hgr 64"));
} // namespace dag
} // namespace kahypar