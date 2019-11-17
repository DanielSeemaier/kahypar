#include "gmock/gmock.h"

#include "dag.h"

#include "kahypar/partition/refinement/acyclic_2way_km1_refiner.h"
#include "kahypar/partition/refinement/acyclic_local_search_refiner.h"
#include "kahypar/partition/refinement/policies/fm_stop_policy.h"

using ::testing::TestWithParam;
using ::testing::Values;
using ::testing::Le;
using ::testing::Eq;

namespace kahypar {
namespace dag {
using Refiner = AcyclicTwoWayKMinusOneRefiner;

class AcyclicTwoWayKM1RefinerTest : public BaseDAGTest, public TestWithParam<const char*> {
 protected:
  void SetUp() override {
    std::string filename;
    PartitionID k;
    std::tie(filename, k) = parseParams();
    loadGraph(filename, k);
  }

  void loadGraph(const std::string& filename, const PartitionID k) override {
    BaseDAGTest::loadGraph(filename, k);
    qg = std::make_unique<AdjacencyMatrixQuotientGraph<DFSCycleDetector>>(hg, context);
    gain_manager = std::make_unique<KMinusOneGainManager>(hg, context);
    refiner = std::make_unique<Refiner>(hg, context, *qg, *gain_manager);
  }

  void runRefiner(const double target_epsilon, std::vector<HypernodeID> refinement_nodes = {}) {
    context.partition.epsilon = target_epsilon;
    context.setupPartWeights(hg.totalWeight());
    if (refinement_nodes.empty()) {
      refinement_nodes = borderNodes();
    }
    const auto initial_km1 = metrics::km1(hg);
    Metrics metrics{metrics::hyperedgeCut(hg), initial_km1, metrics::imbalance(hg, context)};
    refiner->refine(refinement_nodes, {0, 0}, {}, metrics);

    ASSERT_TRUE(qg->isAcyclic());
    ASSERT_TRUE(AdjacencyMatrixQuotientGraph<DFSCycleDetector>(hg, context).isAcyclic());
    ASSERT_THAT(metrics::km1(hg), Le(initial_km1));
    ASSERT_THAT(metrics::km1(hg), Eq(metrics.km1));
    ASSERT_THAT(metrics::imbalance(hg, context), Eq(metrics.imbalance));
    ASSERT_THAT(metrics::imbalance(hg, context), Le(target_epsilon));
    ASSERT_FALSE(gain_manager->hasDelta());
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
  std::unique_ptr<AdjacencyMatrixQuotientGraph<DFSCycleDetector>> qg{};
  std::unique_ptr<KMinusOneGainManager> gain_manager{};
  std::unique_ptr<Refiner> refiner{};
};

TEST_P(AcyclicTwoWayKM1RefinerTest, CanRefineAllBorderNodes) {
  gain_manager->initialize();
  qg->rebuild();
  refiner->initialize(0);
  runRefiner(0.03, borderNodes());
}

TEST_P(AcyclicTwoWayKM1RefinerTest, OnlyUsesBorderNodes) {
  // partition using all nodes
  gain_manager->initialize();
  qg->rebuild();
  refiner->initialize(0);
  runRefiner(0.03, allNodes());

  const HyperedgeWeight expected_km1 = metrics::km1(hg);
  const double expected_imbalance = metrics::imbalance(hg, context);

  // partition using only border nodes
  partitionUsingTopologicalOrdering(context.partition.k);
  gain_manager->initialize();
  qg->rebuild();
  refiner->initialize(0);
  runRefiner(0.03, borderNodes());

  ASSERT_THAT(metrics::km1(hg), Eq(expected_km1));
  ASSERT_THAT(metrics::imbalance(hg, context), Eq(expected_imbalance));
}

TEST_P(AcyclicTwoWayKM1RefinerTest, AnotherTestCase) {
  const auto my_variable = 10;
  if (my_variable >= 5) {
    std::cout << "Hello World" << std::endl;
  }
  const auto f = []() {
    return 0;
  };
  std::cout << f() << std::endl;
  std::cout << "This is a very clean and readable font that every programmer should use for day to day business" << std::endl;

}

TEST_P(AcyclicTwoWayKM1RefinerTest, CanRefineDuringUncoarsening) {
  auto contractions = contractArbitrarily(3, 250, true);

  qg->rebuild();
  gain_manager->initialize();
  refiner->initialize(0);

  for (auto rit = contractions.crbegin(); rit != contractions.crend(); ++rit) {
    refiner->preUncontraction(rit->u);
    qg->preUncontraction(rit->u);
    gain_manager->preUncontraction(rit->u);
    hg.uncontract(*rit);
    gain_manager->postUncontraction(rit->u, {rit->v});
    qg->postUncontraction(rit->u, {rit->v});
    refiner->postUncontraction(rit->u, {rit->v});

    runRefiner(0.03, {rit->u, rit->v});
  }
}

INSTANTIATE_TEST_CASE_P(GRAPH_C17_K_2, AcyclicTwoWayKM1RefinerTest, Values("test_instances/c17.hgr 2"));

INSTANTIATE_TEST_CASE_P(GRAPH_C3540_K_2, AcyclicTwoWayKM1RefinerTest, Values("test_instances/c3540.hgr 2"));
} // namespace dag
} // namespace kahypar