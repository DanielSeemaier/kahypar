#include "gmock/gmock.h"

#include "dag.h"

#include "kahypar/partition/refinement/acyclic_soft_rebalance_refiner.h"

using ::testing::TestWithParam;
using ::testing::Values;
using ::testing::Le;

namespace kahypar {
namespace dag {
class SoftRebalanceTest : public BaseDAGTest, public TestWithParam<const char*> {
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
    refiner = std::make_unique<AcyclicSoftRebalanceRefiner>(hg, context, *qg, *gain_manager);
  }

  void runRefiner(const double target_epsilon) {
    context.partition.epsilon = target_epsilon;
    context.setupPartWeights(hg.totalWeight());
    std::vector<HypernodeID> refinement_nodes{};
    Metrics metrics{metrics::hyperedgeCut(hg), metrics::km1(hg), metrics::imbalance(hg, context)};
    refiner->refine(refinement_nodes, {0, 0}, {}, metrics);

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
  std::unique_ptr<AdjacencyMatrixQuotientGraph<DFSCycleDetector>> qg;
  std::unique_ptr<KMinusOneGainManager> gain_manager;
  std::unique_ptr<AcyclicSoftRebalanceRefiner> refiner;
};

TEST_P(SoftRebalanceTest, CanRebalanceWithoutCrash) {
  applyImbalancedPartition(std::max(context.partition.k / 4, 1), 0.75);
  qg->rebuild();
  gain_manager->initialize();
  refiner->initialize(0);
  runRefiner(0.03);
//  ASSERT_THAT(metrics::imbalance(hg, context), Le(0.03));
}

//TEST_P(HardRebalanceTest, CanRebalanceDuringUncoarsening) {
//  applyImbalancedPartition(std::max(context.partition.k / 4, 1), 0.75);
//  auto contractions = contractArbitrarily(3, 250, true);
//
//  qg->rebuild(); // quotient graph cannot handle contractions
//  gain_manager->initialize(); // gain manager cannot handle contractions
//  refiner->initialize(0);
//
//  const double FINAL_EPSILON = 0.03;
//  const double initial_epsilon = metrics::imbalance(hg, context);
//
//  HypernodeID least_num_nodes = hg.currentNumNodes();
//
//  for (auto rit = contractions.crbegin(); rit != contractions.crend(); ++rit) {
//    refiner->preUncontraction(rit->u);
//    qg->preUncontraction(rit->u);
//    gain_manager->preUncontraction(rit->u);
//    hg.uncontract(*rit);
//    gain_manager->postUncontraction(rit->u, {rit->v});
//    qg->postUncontraction(rit->u, {rit->v});
//    refiner->postUncontraction(rit->u, {rit->v});
//
//    const double progress = 1.0 - (hg.initialNumNodes() - hg.currentNumNodes())
//                                  / static_cast<double>(hg.initialNumNodes() - least_num_nodes);
//    const double delta = initial_epsilon - FINAL_EPSILON;
//    const double epsilon = initial_epsilon - delta * progress;
//
//    runRefiner(epsilon);
//    qg->resetChanged();
//    ASSERT_THAT(metrics::imbalance(hg, context), Le(epsilon));
//  }
//
//  ASSERT_THAT(metrics::imbalance(hg, context), Le(FINAL_EPSILON));
//}

INSTANTIATE_TEST_CASE_P(GRAPH_C17_K_2, SoftRebalanceTest, Values("test_instances/c17.hgr 2"));

INSTANTIATE_TEST_CASE_P(GRAPH_C17_K_4, SoftRebalanceTest, Values("test_instances/c17.hgr 4"));

INSTANTIATE_TEST_CASE_P(GRAPH_C3540_K_2, SoftRebalanceTest, Values("test_instances/c3540.hgr 2"));

INSTANTIATE_TEST_CASE_P(GRAPH_C3540_K_4, SoftRebalanceTest, Values("test_instances/c3540.hgr 4"));

//INSTANTIATE_TEST_CASE_P(GRAPH_C3540_K_32, HardRebalanceTest, Values("test_instances/c3540.hgr 32"));
//
//INSTANTIATE_TEST_CASE_P(GRAPH_C3540_K_64, HardRebalanceTest, Values("test_instances/c3540.hgr 64"));
} // namespace dag
} // namespace kahypar