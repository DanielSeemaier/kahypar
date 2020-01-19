#include "gmock/gmock.h"

#include "dag.h"

#include "kahypar/io/partitioning_output.h"
#include "kahypar/partition/refinement/acyclic_2way_km1_refiner.h"
#include "kahypar/partition/refinement/wip_refiner.h"
#include "kahypar/partition/refinement/acyclic_local_search_refiner.h"
#include "kahypar/partition/refinement/policies/fm_stop_policy.h"
#include "kahypar/definitions.h"
#include "kahypar/io/hypergraph_io.h"
#include "kahypar/partition/coarsening/i_coarsener.h"
#include "kahypar/partition/context.h"
#include "kahypar/partition/factories.h"
#include "kahypar/partition/fixed_vertices.h"
#include "kahypar/partition/multilevel.h"
#include "kahypar/partition/preprocessing/louvain.h"
#include "kahypar/partition/refinement/i_refiner.h"

using namespace kahypar;

TEST(Debug, Debug) {
  Hypergraph hg = dag::loadHypergraph("/Users/danielseemaier/tmp/current_hypergraph.hgr");
  LOG << hg.initialNumNodes();
  std::vector<PartitionID> partition;
  io::readPartitionFile("/Users/danielseemaier/tmp/current_partition.part", partition);
  LOG << partition.size();
  hg.resetPartitioning();
  hg.changeK(2);
  hg.setPartition(partition);

  Context ctx;
  ctx.partition.k = 2;
  ctx.partition.objective = Objective::km1;
  ctx.partition.epsilon = 0.03;
  ctx.setupPartWeights(hg.totalWeight());
  ctx.local_search.fm.max_number_of_fruitless_moves = 350;

  LOG << AdjacencyMatrixQuotientGraph<DFSCycleDetector>(hg, ctx).isAcyclic();
  LOG << metrics::km1(hg);

  AcyclicTwoWayKMinusOneRefiner refiner(hg, ctx);

  AcyclicClusterCoarseningV2<HeavyEdgeScore,
      NoWeightPenalty,
      IgnoreCommunityStructure,
      NormalPartitionPolicy,
      BestRatingPreferringUnmatched<>,
      AllowFreeOnFixedFreeOnFreeFixedOnFixed, double> coarsener(
      hg, ctx, 0);

  coarsener.coarsen(320);
  hg.initializeNumCutHyperedges();
  coarsener.uncoarsen(refiner);

  LOG << metrics::km1(hg);
  LOG << AdjacencyMatrixQuotientGraph<DFSCycleDetector>(hg, ctx).isAcyclic();

}