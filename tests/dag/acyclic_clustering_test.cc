#include "gmock/gmock.h"

#include "dag.h"

#include "kahypar/partition/refinement/acyclic_local_search_refiner.h"
#include "kahypar/partition/refinement/policies/fm_stop_policy.h"
#include "kahypar/dag/acyclic_clustering.h"

using ::testing::TestWithParam;
using ::testing::Values;
using ::testing::Le;
using ::testing::Eq;

namespace kahypar {
namespace dag {
using LocalSearchRefiner = AcyclicLocalSearchRefiner<NumberOfFruitlessMovesStopsSearch>;

class AcyclicClusteringTest : public BaseDAGTest, public TestWithParam<const char*> {
 protected:
  void SetUp() override {
    loadGraph(GetParam(), Hypergraph::kInvalidPartition);
  }

  void loadGraph(const std::string& filename, const PartitionID k) override {
    hg = loadHypergraph(filename);
    for (const HypernodeID& hn : hg.nodes()) {
      hg.setNodePart(hn, 0);
    }
  }

  void contractClustering(const std::vector<PartitionID>& clustering) {
    for (const HypernodeID& hn : hg.nodes()) {
      if (clustering[hn] != hn) {
        hg.contract(clustering[hn], hn);
      }
    }
  }

  void ASSERT_THAT_CLUSTERING_IS_CONTRACTIBLE(const std::vector<PartitionID>& clustering) {
    std::vector<Hypergraph::Memento> contractions;

    for (const HypernodeID& hn : hg.nodes()) {
      if (clustering[hn] != hn) {
        contractions.push_back(hg.contract(clustering[hn], hn));
        ASSERT_TRUE(isAcyclic(hg)) << hg.initialNumNodes() << " / " << contractions.size();
      }
    }

    for (auto rit = contractions.crbegin(); rit != contractions.crend(); ++rit) {
      hg.uncontract(*rit);
    }
  }
};

TEST_P(AcyclicClusteringTest, SingleIterationOfClusteringIsContractible) {
  const auto clustering = findAcyclicClustering(hg, context, 1.0);
  ASSERT_THAT_CLUSTERING_IS_CONTRACTIBLE(clustering);
}

TEST_P(AcyclicClusteringTest, TwoIterationsOfClusteringAreContractible) {
  const auto clustering_iter1 = findAcyclicClustering(hg, context, 1.0);
  contractClustering(clustering_iter1);
  const auto clustering_iter2 = findAcyclicClustering(hg, context, 1.0);
  ASSERT_THAT_CLUSTERING_IS_CONTRACTIBLE(clustering_iter2);
}

TEST_P(AcyclicClusteringTest, MaxClusterWeightIsRespected) {
  const double max_weight_fraction = 0.01;
  const auto max_weight = static_cast<HypernodeWeight>(max_weight_fraction * hg.totalWeight()) + 1;
  HypernodeID last_num_nodes = 0;
  do {
    last_num_nodes = hg.currentNumNodes();
    const auto clustering = findAcyclicClustering(hg, context, max_weight_fraction);
    contractClustering(clustering);
  } while (last_num_nodes > hg.currentNumNodes());

  for (const HypernodeID& hn : hg.nodes()) {
    ASSERT_THAT(hg.nodeWeight(hn), Le(max_weight));
  }
}

//INSTANTIATE_TEST_CASE_P(GRAPH_STAR, AcyclicClusteringTest, Values("test_instances/star.hgr"));

INSTANTIATE_TEST_CASE_P(GRAPH_C17, AcyclicClusteringTest, Values("test_instances/c17.hgr"));

INSTANTIATE_TEST_CASE_P(GRAPH_C880, AcyclicClusteringTest, Values("test_instances/c880.hgr"));

INSTANTIATE_TEST_CASE_P(GRAPH_C7552, AcyclicClusteringTest, Values("test_instances/c7552.hgr"));

//INSTANTIATE_TEST_CASE_P(GRAPH_VIBROBOX, AcyclicClusteringTest, Values("test_instances/vibrobox.hgr"));

//INSTANTIATE_TEST_CASE_P(GRAPH_2MM, AcyclicClusteringTest, Values("test_instances/2mm.hgr"));
} // namespace dag
} // namespace kahypar