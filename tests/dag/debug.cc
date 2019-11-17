#include "gmock/gmock.h"

#include "dag.h"

#include "kahypar/partition/refinement/acyclic_2way_km1_refiner.h"
#include "kahypar/partition/refinement/wip_refiner.h"
#include "kahypar/partition/refinement/acyclic_local_search_refiner.h"
#include "kahypar/partition/refinement/policies/fm_stop_policy.h"

using namespace kahypar;

TEST(Debug, Debug) {
  Hypergraph hg = dag::loadHypergraph("/Users/danielseemaier/HyperDAG/GraphHypergraphs/ISPD98/ibm01.graph.hgr");
  std::vector<PartitionID> partition;
  io::readPartitionFile("/Users/danielseemaier/ibm01.part.2", partition);
  hg.setPartition(partition);
  hg.initializeNumCutHyperedges();

  Context context;
  context.partition.k = 2;
  context.partition.epsilon = 0.03;
  context.setupPartWeights(hg.totalWeight());
  context.local_search.fm.max_number_of_fruitless_moves = 350;




  {
    LOG << "WipRefiner";
    AdjacencyMatrixQuotientGraph<DFSCycleDetector> qg(hg, context);
    KMinusOneGainManager gain_manager(hg, context);
    gain_manager.initialize();
    UncontractionGainChanges changes{}; // dummy
    Metrics current_metrics = {metrics::hyperedgeCut(hg),
                               metrics::km1(hg),
                               metrics::imbalance(hg, context)};
    WipRefiner local_search_refiner(hg, context, qg, gain_manager);
    local_search_refiner.initialize(0);
    std::vector<HypernodeID> local_search_refinement_nodes;
    for (const HypernodeID& hn : hg.nodes()) {
      local_search_refinement_nodes.push_back(hn);
    }
    local_search_refiner.refine(local_search_refinement_nodes, {0, 0}, changes, current_metrics);
    local_search_refiner.printSummary();
  }

  hg.resetPartitioning();
  hg.setPartition(partition);
  hg.initializeNumCutHyperedges();

  {
    LOG << "AcyclicTwoWayKMinusOneRefiner";
    AdjacencyMatrixQuotientGraph<DFSCycleDetector> qg(hg, context);
    KMinusOneGainManager gain_manager(hg, context);
    gain_manager.initialize();
    UncontractionGainChanges changes{}; // dummy
    Metrics current_metrics = {metrics::hyperedgeCut(hg),
                               metrics::km1(hg),
                               metrics::imbalance(hg, context)};
    AcyclicTwoWayKMinusOneRefiner local_search_refiner(hg, context, qg, gain_manager);
    local_search_refiner.initialize(0);
    std::vector<HypernodeID> local_search_refinement_nodes;
    for (const HypernodeID& hn : hg.nodes()) {
      local_search_refinement_nodes.push_back(hn);
    }
    local_search_refiner.refine(local_search_refinement_nodes, {0, 0}, changes, current_metrics);
    local_search_refiner.printSummary();
  }
}