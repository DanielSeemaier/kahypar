#include "gmock/gmock.h"

#include "dag.h"

#include "kahypar/partition/refinement/acyclic_2way_km1_refiner.h"
#include "kahypar/partition/refinement/wip_refiner.h"
#include "kahypar/partition/refinement/acyclic_local_search_refiner.h"
#include "kahypar/partition/refinement/policies/fm_stop_policy.h"

using namespace kahypar;

TEST(Debug, Debug) {
  Hypergraph hg = dag::loadHypergraph("/Users/danielseemaier/tmp/syr2k.hgr");
  LOG << hg.initialNumNodes();
  std::vector<PartitionID> partition;
  io::readPartitionFile("/Users/danielseemaier/tmp/original_partition.part", partition);
  LOG << partition.size();
  hg.resetPartitioning();
  hg.changeK(32);
  hg.setPartition(partition);

  Context ctx;
  ctx.partition.k = 32;
  ctx.partition.epsilon = 0.03;
  ctx.setupPartWeights(hg.totalWeight());

  LOG << AdjacencyMatrixQuotientGraph<DFSCycleDetector>(hg, ctx).isAcyclic();

}