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

static void _partitionUsingTopologicalOrdering(Hypergraph& hg, const PartitionID k) {
  auto ordering = calculateTopologicalOrdering(hg);
  PartitionID part = 0;
  HypernodeID nodes_per_part = hg.currentNumPins() / k;
  HypernodeID nodes_in_cur_part = 0;

  hg.changeK(k);
  for (const HypernodeID& hn : ordering) {
    hg.setNodePart(hn, part);
    if (nodes_in_cur_part == nodes_per_part) {
      ++part;
      nodes_in_cur_part = 0;
    }
  }
}

static Context _createContext(PartitionID k) {
  Context ctx;
  ctx.partition.k = k;
  return ctx;
}

template<typename Detector>
static void _performMonkeyMoves(Hypergraph& hg, QuotientGraph<Detector>& qg, const std::size_t k, const std::size_t N) {
  Randomize& random = Randomize::instance();
  for (std::size_t i = 0; i < N; ++i) {
    const HypernodeID hn = random.getRandomInt(0, hg.currentNumNodes() - 1);
    const PartitionID from = hg.partID(hn);
    const PartitionID to = random.getRandomInt(0, k - 1);
    if (from != to) { // otherwise changeNodePart crashes
      if (qg.update(hn, from, to)) { // expect behaviour: shouldn't crash ...
        hg.changeNodePart(hn, from, to);
      }
    }
  }
}

TEST(QUOTIENT_GRAPH, HasCorrectSize) {
  Hypergraph hg = _loadHypergraph("test_instances/acyclic.hgr");
  _applyAcyclic4Partition(hg);
  QuotientGraph<DFSCycleDetector> qg(hg, _createContext(4));

  ASSERT_EQ(qg.numberOfNodes(), 4);
  ASSERT_EQ(qg.numberOfEdges(), 3);
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

TEST(QUOTIENT_GRAPH, MixedLegalAndIllegalMoves) {
  Hypergraph hg = _loadHypergraph("test_instances/acyclic.hgr");
  _applyAcyclic4Partition(hg);
  QuotientGraph<DFSCycleDetector> qg(hg, _createContext(4));

  ASSERT_FALSE(qg.update(5, 3, 0));
  ASSERT_TRUE(qg.update(5, 3, 2));
  ASSERT_FALSE(qg.update(5, 2, 1));
  ASSERT_TRUE(qg.update(5, 2, 3));
}

TEST(QUOTIENT_GRAPH, NullMovesAreLegal) {
  Hypergraph hg = _loadHypergraph("test_instances/acyclic.hgr");
  _applyAcyclic4Partition(hg);
  QuotientGraph<DFSCycleDetector> qg(hg, _createContext(4));

  for (const HypernodeID& hn : hg.nodes()) {
    ASSERT_TRUE(qg.update(hn, hg.partID(hn), hg.partID(hn)));
  }
}

TEST(QUOTIENT_GRAPH, MonkeyMovesDontCrashAndRemainsAcyclic) {
  Hypergraph hg = _loadHypergraph("test_instances/acyclic.hgr");
  _applyAcyclic4Partition(hg);
  QuotientGraph<DFSCycleDetector> qg(hg, _createContext(4));

  constexpr std::size_t N = 128;
  constexpr PartitionID k = 4;
  _performMonkeyMoves(hg, qg, k, N);
  ASSERT_TRUE(qg.isAcyclic());
}

TEST(QUOTIENT_GRAPH, C880_K128_InitiallyAcyclic) {
  constexpr PartitionID k = 128;
  Hypergraph hg = _loadHypergraph("test_instances/c880.hgr");
  _partitionUsingTopologicalOrdering(hg, k);
  QuotientGraph<DFSCycleDetector> qg(hg, _createContext(k));
  ASSERT_TRUE(qg.isAcyclic());
}

TEST(QUOTIENT_GRAPH, C880_K64_MonkeyMovesDontCrashAndRemainsAcyclic) {
  constexpr PartitionID k = 64;
  constexpr std::size_t N = 1024;
  Hypergraph hg = _loadHypergraph("test_instances/c880.hgr");
  _partitionUsingTopologicalOrdering(hg, k);
  QuotientGraph<DFSCycleDetector> qg(hg, _createContext(k));
  _performMonkeyMoves(hg, qg, k, N);
  ASSERT_TRUE(qg.isAcyclic());
}

TEST(QUOTIENT_GRAPH, C880_K128_MonkeyMovesDontCrashAndRemainsAcyclic) {
  constexpr PartitionID k = 128;
  constexpr std::size_t N = 1024;
  Hypergraph hg = _loadHypergraph("test_instances/c880.hgr");
  _partitionUsingTopologicalOrdering(hg, k);
  QuotientGraph<DFSCycleDetector> qg(hg, _createContext(k));
  _performMonkeyMoves(hg, qg, k, N);
  ASSERT_TRUE(qg.isAcyclic());
}
} // namespace dag
} // namespace kahypar