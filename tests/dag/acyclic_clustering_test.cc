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

class AcyclicClusteringTest : public BaseDAGTest, public TestWithParam<const char *> {
protected:
  void SetUp() override {
    loadGraph(GetParam(), Hypergraph::kInvalidPartition);
  }

  void loadGraph(const std::string &filename, const PartitionID k) override {
    hg = loadHypergraph(filename);
    for (const HypernodeID &hn : hg.nodes()) {
      hg.setNodePart(hn, 0);
    }
  }

  void contractClustering(const std::vector<HypernodeID> &clustering) {
    for (const HypernodeID &hn : hg.nodes()) {
      if (clustering[hn] != hn) {
        hg.contract(clustering[hn], hn);
      }
    }
  }

  void ASSERT_THAT_CLUSTERING_IS_CONTRACTIBLE(const std::vector<HypernodeID> &clustering) {
    std::vector<Hypergraph::Memento> contractions;

    for (const HypernodeID &hn : hg.nodes()) {
      if (clustering[hn] != hn) {
        contractions.push_back(hg.contract(clustering[hn], hn));
//        LOG << "contractions.push_back(hg.contract(clustering[" << hn << "], " << hn << "));";
        const bool acyclic = isAcyclic(hg);
        if (!acyclic) {
          const auto cycle = _debug_findCycleContainingHN(hg, clustering[hn]);
          LOG << "Acyclic:" << cycle;
          std::exit(1);
        }
        ASSERT_TRUE(acyclic);
      }
    }

    for (auto rit = contractions.crbegin(); rit != contractions.crend(); ++rit) {
      hg.uncontract(*rit);
    }
  }

  void ASSERT_THAT_CLUSTERING_IS_CONTRACTIBLE2(const std::vector<HypernodeID> &clustering) {
    std::vector<Hypergraph::Memento> contractions;
    const auto top = dag::calculateToplevelValues(hg);

//    hg.printDirectedHypergraphNode(31583);
//    hg.printDirectedHypergraphNode(1944);
//    hg.printDirectedHypergraphNode(31576);
//    hg.printDirectedHypergraphNode(1969);
//    hg.printDirectedHypergraphNode(31602);
//    hg.printDirectedHypergraphNode(31606);

    LOG << V(top[31583]);
    LOG << V(top[1944]);
    LOG << V(top[31576]);
    LOG << V(top[1969]);
    LOG << V(top[31602]);
    LOG << V(top[31606]);
//
//    LOG << clustering[434] << clustering[443];
//    LOG << top[434] << top[clustering[434]] << top[443] << top[clustering[443]];

    LOG << "Contract" << clustering[1944] << 1944;
    contractions.push_back(hg.contract(clustering[1944], 1944));

    LOG << "Contract" << clustering[1969] << 1969;
    contractions.push_back(hg.contract(clustering[1969], 1969));

    LOG << "Contract" << clustering[31602] << 31602;
    LOG << "Contract" << clustering[31606] << 31606;
    contractions.push_back(hg.contract(clustering[31602], 31602));
    contractions.push_back(hg.contract(clustering[31606], 31606));


    const bool acyclic = isAcyclic(hg);
    if (!acyclic) {
      const auto cycle = _debug_findCycleContainingHN(hg, 31583);
      LOG << "Acyclic:" << cycle;
    }
    ASSERT_TRUE(acyclic);

  }

  std::vector<HypernodeID> _debug_findCycleContainingHN(const Hypergraph &hg, const HypernodeID hn) const {
    std::vector<bool> marked(hg.initialNumNodes());
    marked[hn] = true;
    std::vector<HypernodeID> cycle;
    if (_debug_findCycleContainingHN(hg, hn, hn, marked, cycle)) {
      cycle.push_back(hn);
    }
    return cycle;
  }

  bool _debug_findCycleContainingHN(const Hypergraph &hg, const HypernodeID hn, const HypernodeID start,
                                    std::vector<bool> &marked, std::vector<HypernodeID> &cycle) const {
    for (const HyperedgeID &he : hg.incidentTailEdges(hn)) {
      for (const HyperedgeID &head : hg.heads(he)) {
        if (head == start) {
          cycle.push_back(start);
          return true;
        }
        if (marked[head]) {
          continue;
        }

        marked[head] = true;
        if (_debug_findCycleContainingHN(hg, head, start, marked, cycle)) {
          cycle.push_back(head);
        }
      }
    }

    return false;
  }
};

//TEST_P(AcyclicClusteringTest, SingleIterationOfClusteringIsContractible) {
//  const auto clustering = findAcyclicClustering(hg, context, 1.0);
//  ASSERT_THAT_CLUSTERING_IS_CONTRACTIBLE(clustering);
//}

//TEST_P(AcyclicClusteringTest, SingleIterationOfClusteringWithCCIsContractible) {
//  const auto clustering = findAcyclicClusteringWithCycleDetection(hg, context, 1.0);
//  ASSERT_THAT_CLUSTERING_IS_CONTRACTIBLE(clustering);
//}

//TEST_P(AcyclicClusteringTest, TwoIterationsOfClusteringAreContractible) {
//  const auto clustering_iter1 = findAcyclicClustering(hg, context, 1.0);
//  contractClustering(clustering_iter1);
//  const auto clustering_iter2 = findAcyclicClustering(hg, context, 1.0);
//  ASSERT_THAT_CLUSTERING_IS_CONTRACTIBLE(clustering_iter2);
//}

TEST_P(AcyclicClusteringTest, ThreeIterationsOfClusteringWithCCAreContractible) {
  const auto clustering_iter1 = findAcyclicClusteringWithCycleDetection(hg, context, 1.0);
  ASSERT_THAT_CLUSTERING_IS_CONTRACTIBLE(clustering_iter1);
  contractClustering(clustering_iter1);
  LOG << "Level 1 OK";

  const auto clustering_iter2 = findAcyclicClusteringWithCycleDetection(hg, context, 1.0);
  ASSERT_THAT_CLUSTERING_IS_CONTRACTIBLE(clustering_iter2);
  contractClustering(clustering_iter2);
  LOG << "Level 2 OK";

  const auto clustering_iter3 = findAcyclicClusteringWithCycleDetection(hg, context, 1.0);
  ASSERT_THAT_CLUSTERING_IS_CONTRACTIBLE(clustering_iter3);
  LOG << "Level 3 OK";
}

//TEST_P(AcyclicClusteringTest, Debug__RemoveMe) {
//  std::vector<PartitionID> ip;
//  io::readPartitionFile("/Users/danielseemaier/Projects/kahypar/tmp/ibm04.ip", ip);
//  hg.setPartition(ip);
//
//  auto clustering_iter1 = findAcyclicClusteringWithCycleDetection(hg, context, 1.0);
////  for (HypernodeID hn = 0; hn < 1940; ++hn) {
////    clustering_iter1[hn] = hn;
////  }
////  for (HypernodeID hn = 1970; hn < 31600; ++hn) {
////    clustering_iter1[hn] = hn;
////  }
//  ASSERT_THAT_CLUSTERING_IS_CONTRACTIBLE(clustering_iter1);
//}

//TEST_P(AcyclicClusteringTest, MaxClusterWeightIsRespected) {
//  const double max_weight_fraction = 0.01;
//  const auto max_weight = static_cast<HypernodeWeight>(max_weight_fraction * hg.totalWeight()) + 1;
//  HypernodeID last_num_nodes = 0;
//  do {
//    last_num_nodes = hg.currentNumNodes();
//    const auto clustering = findAcyclicClustering(hg, context, max_weight_fraction);
//    contractClustering(clustering);
//  } while (last_num_nodes > hg.currentNumNodes());
//
//  for (const HypernodeID& hn : hg.nodes()) {
//    ASSERT_THAT(hg.nodeWeight(hn), Le(max_weight));
//  }
//}

//INSTANTIATE_TEST_CASE_P(GRAPH_STAR, AcyclicClusteringTest, Values("test_instances/star.hgr"));

//INSTANTIATE_TEST_CASE_P(GRAPH_C17, AcyclicClusteringTest, Values("test_instances/c17.hgr"));

//INSTANTIATE_TEST_CASE_P(GRAPH_C880, AcyclicClusteringTest, Values("test_instances/c880.hgr"));
//
//INSTANTIATE_TEST_CASE_P(GRAPH_C7552, AcyclicClusteringTest, Values("test_instances/c7552.hgr"));
//
INSTANTIATE_TEST_CASE_P(GRAPH_VIBROBOX, AcyclicClusteringTest, Values("test_instances/vibrobox.hgr"));
INSTANTIATE_TEST_CASE_P(GRAPH_2MM, AcyclicClusteringTest, Values("test_instances/2mm_graph.hgr"));
INSTANTIATE_TEST_CASE_P(HYPERGRAPH_IBM01, AcyclicClusteringTest, Values("test_instances/ibm01_hypergraph.hgr"));
INSTANTIATE_TEST_CASE_P(HYPERGRAPH_IBM04, AcyclicClusteringTest, Values("test_instances/ibm04_hypergraph.hgr"));
} // namespace dag
} // namespace kahypar