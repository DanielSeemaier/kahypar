#include "gmock/gmock.h"

#include "kahypar/dag/quotient_graph.h"

#include "dag.h"

using testing::Test;
using testing::Eq;
using testing::Le;
using testing::Lt;

namespace kahypar {
namespace dag {
class QuotientGraphTest : public Test {
 protected:
  void SetUp() override {
    hg = loadHypergraph("test_instances/acyclic.hgr");
  }

  void applyAcyclic4Partition() {
    ctx.partition.k = 4;
    hg.changeK(4);
    hg.setNodePart(6, 0);
    hg.setNodePart(0, 0);
    hg.setNodePart(1, 1);
    hg.setNodePart(2, 2);
    hg.setNodePart(3, 2);
    hg.setNodePart(4, 2);
    hg.setNodePart(5, 3);
  }

  void applySingletonPartition() {
    ctx.partition.k = hg.currentNumNodes();
    hg.changeK(hg.currentNumNodes());
    for (const HypernodeID& hn : hg.nodes()) {
      hg.setNodePart(hn, hn);
    }
  }

  void partitionUsingTopologicalOrdering(const PartitionID k) {
    ctx.partition.k = k;

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

  QuotientGraph<DFSCycleDetector> createQuotientGraph() const {
    return QuotientGraph<DFSCycleDetector>(hg, ctx);
  }

  template<typename Detector>
  void performMonkeyMoves(QuotientGraph<Detector>& qg, const std::size_t k, const std::size_t N) {
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

  Hypergraph hg;
  Context ctx;
};

TEST_F(QuotientGraphTest, HasCorrectSize) {
  applyAcyclic4Partition();
  const auto qg = createQuotientGraph();
  ASSERT_EQ(qg.numberOfNodes(), 4);
  ASSERT_EQ(qg.numberOfEdges(), 3);
}

TEST_F(QuotientGraphTest, AcyclicInitialPartitionHasAcyclicQuotientGraph) {
  applyAcyclic4Partition();
  const auto qg = createQuotientGraph();
  ASSERT_TRUE(qg.isAcyclic());
}

TEST_F(QuotientGraphTest, SimpleLegalMovesAreLegal) {
  applyAcyclic4Partition();
  auto qg = createQuotientGraph();
  ASSERT_TRUE(qg.update(2, 2, 1)); // move 2 from part 2 to part 1
  hg.changeNodePart(2, 2, 1);
  ASSERT_TRUE(qg.update(2, 1, 2)); // revert movement
  hg.changeNodePart(2, 1, 2);
}

TEST_F(QuotientGraphTest, SimpleIllegalMovesAreIllegal) {
  applyAcyclic4Partition();
  auto qg = createQuotientGraph();
  ASSERT_FALSE(qg.update(1, 1, 3)); // move 1 from part 1 to part 3
  ASSERT_FALSE(qg.update(6, 0, 3)); // move 6 from part 0 to part 3
}

TEST_F(QuotientGraphTest, MixedLegalAndIllegalMoves) {
  applyAcyclic4Partition();
  auto qg = createQuotientGraph();
  ASSERT_FALSE(qg.update(5, 3, 0));
  ASSERT_TRUE(qg.update(5, 3, 2));
  ASSERT_FALSE(qg.update(5, 2, 1));
  ASSERT_TRUE(qg.update(5, 2, 3));
}

TEST_F(QuotientGraphTest, NullMovesAreLegal) {
  applyAcyclic4Partition();
  auto qg = createQuotientGraph();
  for (const HypernodeID& hn : hg.nodes()) {
    ASSERT_TRUE(qg.update(hn, hg.partID(hn), hg.partID(hn)));
  }
}

TEST_F(QuotientGraphTest, MonkeyMovesDontCrashAndRemainsAcyclic) {
  applyAcyclic4Partition();
  auto qg = createQuotientGraph();
  constexpr std::size_t N = 128;
  constexpr PartitionID k = 4;
  performMonkeyMoves(qg, k, N);
  ASSERT_TRUE(qg.isAcyclic());
}

TEST_F(QuotientGraphTest, C880_K128_InitiallyAcyclic) {
  constexpr PartitionID k = 128;
  hg = loadHypergraph("test_instances/c880.hgr");
  partitionUsingTopologicalOrdering(k);
  const auto qg = createQuotientGraph();
  ASSERT_TRUE(qg.isAcyclic());
}

TEST_F(QuotientGraphTest, C880_K64_MonkeyMovesDontCrashAndRemainsAcyclic) {
  constexpr PartitionID k = 64;
  constexpr std::size_t N = 1024;
  hg = loadHypergraph("test_instances/c880.hgr");
  partitionUsingTopologicalOrdering(k);
  auto qg = createQuotientGraph();
  performMonkeyMoves(qg, k, N);
  ASSERT_TRUE(qg.isAcyclic());
}

TEST_F(QuotientGraphTest, C880_K128_MonkeyMovesDontCrashAndRemainsAcyclic) {
  constexpr PartitionID k = 128;
  constexpr std::size_t N = 1024;
  hg = loadHypergraph("test_instances/c880.hgr");
  partitionUsingTopologicalOrdering(k);
  auto qg = createQuotientGraph();
  performMonkeyMoves(qg, k, N);
  ASSERT_TRUE(qg.isAcyclic());
}

TEST_F(QuotientGraphTest, C17_WeakTopologicalOrderingIsCorrect) {
  hg = loadHypergraph("test_instances/c17.hgr");
  applySingletonPartition();
  auto qg = createQuotientGraph();
  auto ordering = qg.computeWeakTopologicalOrdering();
  ASSERT_THAT(ordering[0], Eq(2));
  ASSERT_THAT(ordering[1], Eq(2));
  ASSERT_THAT(ordering[2], Eq(1));
  ASSERT_THAT(ordering[3], Eq(3));
  ASSERT_THAT(ordering[4], Eq(0));
  ASSERT_THAT(ordering[5], Eq(1));
  ASSERT_THAT(ordering[6], Eq(0));
  ASSERT_THAT(ordering[7], Eq(0));
  ASSERT_THAT(ordering[8], Eq(0));
  ASSERT_THAT(ordering[9], Eq(3));
  ASSERT_THAT(ordering[10], Eq(0));
}

TEST_F(QuotientGraphTest, C3540_StrictTopologicalOrderingHolds) {
  hg = loadHypergraph("test_instances/c3540.hgr");
  applySingletonPartition();
  auto qg = createQuotientGraph();
  auto strict = qg.computeTopologicalOrdering();

  for (const HyperedgeID& he : hg.edges()) {
    for (const HypernodeID& hh : hg.heads(he)) {
      for (const HypernodeID& ht : hg.tails(he)) {
        ASSERT_THAT(strict[ht], Lt(strict[hh]));
      }
    }
  }
}

TEST_F(QuotientGraphTest, C3540_WeakTopologicalOrderingIsLessThanStrictTopologicalOrdering) {
  hg = loadHypergraph("test_instances/c3540.hgr");
  applySingletonPartition();
  auto qg = createQuotientGraph();
  auto weak = qg.computeWeakTopologicalOrdering();
  auto strict = qg.computeTopologicalOrdering();
  for (const HypernodeID& hn : hg.nodes()) {
    ASSERT_THAT(weak[hn], Le(strict[hn]));
  }
}

static void ASSERT_THAT_QMG_MOVES_ARE_LEGAL(Hypergraph& hg, QuotientGraph<DFSCycleDetector>& qg) {
  auto qmg = qg.computeMoveGraph();
  for (const HypernodeID& from : hg.nodes()) {
    for (const QNodeID& to : qmg.outs(from)) {
      bool move_from_to = qg.update(from, from, to);
      ASSERT_THAT(move_from_to, Eq(true));
      hg.changeNodePart(from, from, to);
      bool move_to_from = qg.update(from, to, from);
      ASSERT_THAT(move_to_from, Eq(true));
      hg.changeNodePart(from, to, from);
    }
  }
}

TEST_F(QuotientGraphTest, C17_MovesAlongMoveGraphAreLegal) {
  hg = loadHypergraph("test_instances/c17.hgr");
  applySingletonPartition();
  auto qg = createQuotientGraph();
  qg.reduceToOneRoot();
  ASSERT_THAT_QMG_MOVES_ARE_LEGAL(hg, qg);
}

TEST_F(QuotientGraphTest, C3540_MovesAlongMoveGraphAreLegal) {
  hg = loadHypergraph("test_instances/c3540.hgr");
  applySingletonPartition();
  auto qg = createQuotientGraph();
  qg.reduceToOneRoot();
  ASSERT_THAT_QMG_MOVES_ARE_LEGAL(hg, qg);
}

static void ASSERT_THAT_QMG_IS_SCC(const QuotientMoveGraph& qmg) {
  for (QNodeID u = 0; u < qmg.numberOfNodes(); ++u) {
    std::vector<bool> connected_to(qmg.numberOfNodes());
    std::vector<QNodeID> todo{u};

    while (!todo.empty()) {
      QNodeID from = todo.back();
      todo.pop_back();
      if (!connected_to[from]) {
        connected_to[from] = true;
        const auto& outs = qmg.outs(from);
        todo.insert(todo.end(), outs.begin(), outs.end());
      }
    }

    bool connected_to_all = std::all_of(connected_to.begin(), connected_to.end(), [](auto value) { return value; });
    ASSERT_THAT(connected_to_all, Eq(true)) << "bad node: " << u;
  }
}

TEST_F(QuotientGraphTest, C17_MoveGraphIsSCC) {
  hg = loadHypergraph("test_instances/c17.hgr");
  applySingletonPartition();
  auto qg = createQuotientGraph();
  qg.reduceToOneRoot();
  const auto qmg = qg.computeMoveGraph();
  ASSERT_THAT_QMG_IS_SCC(qmg);
}

TEST_F(QuotientGraphTest, C3540_MoveGraphIsSCC) {
  hg = loadHypergraph("test_instances/c3540.hgr");

  { // test with singleton partitions
    applySingletonPartition();
    auto qg = createQuotientGraph();
    qg.reduceToOneRoot();
    const auto qmg = qg.computeMoveGraph();
    ASSERT_THAT_QMG_IS_SCC(qmg);
  }

  for (const PartitionID& k : {2, 4, 8, 16, 32, 64, 128, 256}) {
    hg.resetPartitioning();
    partitionUsingTopologicalOrdering(k);
    auto qg = createQuotientGraph();
    qg.reduceToOneRoot();
    const auto qmg = qg.computeMoveGraph();
    ASSERT_THAT_QMG_IS_SCC(qmg);
  }
}
} // namespace dag
} // namespace kahypar