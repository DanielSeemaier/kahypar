#include <kahypar/partition/refinement/policies/fm_stop_policy.h>
#include "gmock/gmock.h"

#include "dag.h"

#include "kahypar/partition/refinement/acyclic_local_search_refiner.h"

using ::testing::TestWithParam;
using ::testing::Values;
using ::testing::Le;
using ::testing::Eq;

namespace kahypar {
namespace dag {
using LocalSearchRefiner = AcyclicLocalSearchRefiner<NumberOfFruitlessMovesStopsSearch>;

class LocalSearchTest : public BaseDAGTest, public TestWithParam<const char*> {
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
    refiner = std::make_unique<LocalSearchRefiner>(hg, context, *qg, *gain_manager);
  }

  std::vector<HypernodeID> borderNodes() const {
    std::vector<HypernodeID> nodes;
    for (const HypernodeID& hn : hg.nodes()) {
      if (hg.isBorderNode(hn)) {
        nodes.push_back(hn);
      }
    }
    return nodes;
  }

  std::vector<HypernodeID> allNodes() const {
    std::vector<HypernodeID> nodes;
    for (const HypernodeID& hn : hg.nodes()) {
      nodes.push_back(hn);
    }
    return nodes;
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
  std::unique_ptr<LocalSearchRefiner> refiner{};
};

TEST_P(LocalSearchTest, CanRefineAllBorderNodes) {
  gain_manager->initialize();
  qg->rebuild();
  refiner->initialize(0);
  runRefiner(0.03, borderNodes());
}

TEST_P(LocalSearchTest, OnlyUsesBorderNodes) {
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

TEST_P(LocalSearchTest, CanRefineDuringUncoarsening) {
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

INSTANTIATE_TEST_CASE_P(GRAPH_C17_K_2, LocalSearchTest, Values("test_instances/c17.hgr 2"));

INSTANTIATE_TEST_CASE_P(GRAPH_C17_K_4, LocalSearchTest, Values("test_instances/c17.hgr 4"));

INSTANTIATE_TEST_CASE_P(GRAPH_C3540_K_4, LocalSearchTest, Values("test_instances/c3540.hgr 4"));

INSTANTIATE_TEST_CASE_P(GRAPH_C3540_K_32, LocalSearchTest, Values("test_instances/c3540.hgr 32"));

INSTANTIATE_TEST_CASE_P(GRAPH_C3540_K_64, LocalSearchTest, Values("test_instances/c3540.hgr 64"));
} // namespace dag
} // namespace kahypar