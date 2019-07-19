#include "gmock/gmock.h"

#include "dag.h"

#include "kahypar/partition/refinement/acyclic_hard_rebalance_refiner.h"

using ::testing::TestWithParam;
using ::testing::Values;
using ::testing::Le;

namespace kahypar {
namespace dag {
class HardRebalanceTest : public BaseDAGTest, public TestWithParam<const char*> {
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
    refiner = std::make_unique<AcyclicHardRebalanceRefiner>(hg, context, *qg, *gain_manager);
  }

  // assigns the first num_blocks blocks fraction*\V| nodes
  void applyImbalancedPartition(const PartitionID num_blocks, const double fraction) {
    const auto& ordering = calculateTopologicalOrdering(hg);
    hg.resetPartitioning();

    const HypernodeID size_overloaded_block = fraction * hg.initialNumNodes() / num_blocks;
    const HypernodeID size_underloaded_block = (1.0 - fraction) * hg.initialNumNodes() / (context.partition.k - num_blocks);

    PartitionID cur_part = 0;
    HypernodeID cur_size = 0;
    HypernodeID cur_hn = 0;
    while (cur_hn < hg.initialNumNodes()) {
      hg.setNodePart(ordering[cur_hn], cur_part);
      ++cur_size;
      ++cur_hn;

      if (cur_part < num_blocks && cur_size >= size_overloaded_block && cur_part < context.partition.k - 1) {
        cur_size = 0;
        ++cur_part;
      } else if (cur_part >= num_blocks && cur_size >= size_underloaded_block && cur_part < context.partition.k - 1) {
        cur_size = 0;
        ++cur_part;
      }
    }

    qg->rebuild();
    gain_manager->initialize();
    ASSERT(AdjacencyMatrixQuotientGraph<DFSCycleDetector>(hg, context).isAcyclic());
  }

  void runRefiner(const double target_epsilon) {
    context.partition.epsilon = target_epsilon;
    context.setupPartWeights(hg.totalWeight());
    std::vector<HypernodeID> refinement_nodes{};
    Metrics metrics{metrics::hyperedgeCut(hg), metrics::km1(hg), metrics::imbalance(hg, context)};
    refiner->refine(refinement_nodes, {0, 0}, {}, metrics);
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
  std::unique_ptr<AcyclicHardRebalanceRefiner> refiner;
};

TEST_P(HardRebalanceTest, CanRebalanceWithoutCrash) {
  applyImbalancedPartition(std::max(context.partition.k / 4, 1), 0.75);
  refiner->initialize(0);
  runRefiner(0.03);
  ASSERT_THAT(metrics::imbalance(hg, context), Le(0.03));
}

TEST_P(HardRebalanceTest, CanRebalanceDuringUncoarsening) {
  applyImbalancedPartition(std::max(context.partition.k / 4, 1), 0.75);
  auto contractions = contractArbitrarily(3, 250, true);

  qg->rebuild(); // quotient graph cannot handle contractions
  gain_manager->initialize(); // gain manager cannot handle contractions
  refiner->initialize(0);

  const double FINAL_EPSILON = 0.03;
  const double initial_epsilon = metrics::imbalance(hg, context);

  HypernodeID least_num_nodes = hg.currentNumNodes();

  for (auto rit = contractions.crbegin(); rit != contractions.crend(); ++rit) {
    refiner->preUncontraction(rit->u);
    qg->preUncontraction(rit->u);
    gain_manager->preUncontraction(rit->u);
    hg.uncontract(*rit);
    gain_manager->postUncontraction(rit->u, {rit->v});
    qg->postUncontraction(rit->u, {rit->v});
    refiner->postUncontraction(rit->u, {rit->v});

    const double progress = 1.0 - (hg.initialNumNodes() - hg.currentNumNodes())
                                  / static_cast<double>(hg.initialNumNodes() - least_num_nodes);
    const double delta = initial_epsilon - FINAL_EPSILON;
    const double epsilon = initial_epsilon - delta * progress;

    runRefiner(epsilon);
    qg->resetChanged();
    ASSERT_THAT(metrics::imbalance(hg, context), Le(epsilon));
  }

  ASSERT_THAT(metrics::imbalance(hg, context), Le(FINAL_EPSILON));
}

INSTANTIATE_TEST_CASE_P(GRAPH_C17_K_2, HardRebalanceTest, Values("test_instances/c17.hgr 2"));

INSTANTIATE_TEST_CASE_P(GRAPH_C17_K_4, HardRebalanceTest, Values("test_instances/c17.hgr 4"));

INSTANTIATE_TEST_CASE_P(GRAPH_C3540_K_2, HardRebalanceTest, Values("test_instances/c3540.hgr 2"));

INSTANTIATE_TEST_CASE_P(GRAPH_C3540_K_4, HardRebalanceTest, Values("test_instances/c3540.hgr 4"));

//INSTANTIATE_TEST_CASE_P(GRAPH_C3540_K_32, HardRebalanceTest, Values("test_instances/c3540.hgr 32"));
//
//INSTANTIATE_TEST_CASE_P(GRAPH_C3540_K_64, HardRebalanceTest, Values("test_instances/c3540.hgr 64"));
} // namespace dag
} // namespace kahypar