#include "gmock/gmock.h"

#include "kahypar/dag/quotient_graph.h"

#include "dag.h"

using ::testing::Eq;
using ::testing::ContainerEq;

namespace kahypar {
namespace dag {
static void _applyAcyclic4Partition(Hypergraph &hg) {
  hg.changeK(4);
  hg.setNodePart(6, 0);
  hg.setNodePart(0, 0);
  hg.setNodePart(1, 1);
  hg.setNodePart(2, 2);
  hg.setNodePart(3, 2);
  hg.setNodePart(4, 2);
  hg.setNodePart(5, 3);
}

/*static void _applyCyclic4Partition(Hypergraph &hg) {
  hg.changeK(4);
  hg.setNodePart(6, 0);
  hg.setNodePart(0, 0);
  hg.setNodePart(1, 2);
  hg.setNodePart(2, 1);
  hg.setNodePart(3, 2);
  hg.setNodePart(4, 2);
  hg.setNodePart(5, 3);
}*/

static Context _createContext(PartitionID k) {
  Context ctx;
  ctx.partition.k = k;
  return ctx;
}

TEST(QUOTIENT_GRAPH, AcyclicInitialPartitionHasAcyclicQuotientGraph) {
  Hypergraph hg = _loadHypergraph("test_instances/acyclic.hgr");
  _applyAcyclic4Partition(hg);
  QuotientGraph<DFSCycleDetector> qg(hg, _createContext(4));

  ASSERT_TRUE(qg.isAcyclic());
}

TEST(QUOTIENT_GRAPH, SimpleLegalMovesAreLegal) {
  Hypergraph hg = _loadHypergraph("test_instances/acyclic.hgr");
  _applyAcyclic4Partition(hg);
  QuotientGraph<DFSCycleDetector> qg(hg, _createContext(4));

  ASSERT_TRUE(qg.update(2, 2, 1)); // move 2 from part 2 to part 1
  hg.changeNodePart(2, 2, 1);
  ASSERT_TRUE(qg.update(2, 1, 2)); // revert movement
  hg.changeNodePart(2, 1, 2);
}

TEST(QUOTIENT_GRAPH, SimpleIllegalMovesAreIllegal) {
  Hypergraph hg = _loadHypergraph("test_instances/acyclic.hgr");
  _applyAcyclic4Partition(hg);
  QuotientGraph<DFSCycleDetector> qg(hg, _createContext(4));

  ASSERT_FALSE(qg.update(1, 1, 3)); // move 1 from part 1 to part 3
  ASSERT_FALSE(qg.update(6, 0, 3)); // move 6 from part 0 to part 3
}

// not supported (and not needed) for now
/*TEST(QUOTIENT_GRAPH, CyclicInitialPartitionHasCyclicQuotientGraph) {
  Hypergraph  hg = _loadHypergraph("test_instances/acyclic.hgr");
  _applyCyclic4Partition(hg);
  QuotientGraph<DFSCycleDetector> qg(hg, _createContext(4));
  ASSERT_FALSE(qg.isAcyclic());
}*/
} // namespace dag
} // namespace kahypar