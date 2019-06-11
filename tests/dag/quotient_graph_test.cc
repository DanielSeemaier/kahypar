#include "gmock/gmock.h"

#include "kahypar/dag/quotient_graph.h"

#include "dag.h"

using testing::Test;
using testing::Eq;
using testing::Le;
using testing::Lt;
using testing::Gt;

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

    hg.resetPartitioning();
    hg.changeK(k);
    for (const HypernodeID& hn : ordering) {
      hg.setNodePart(hn, part);
      ++nodes_in_cur_part;
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

  void ASSERT_THAT_MOVES_ALONG_MOVE_PATH_ARE_LEGAL(const std::string& filename) {
    hg = loadHypergraph(filename);
    applySingletonPartition();

    auto qg = createQuotientGraph();
    qg.reduceToOneRoot();

    // find partition that has the longest distance to root
    auto ordering = qg.computeWeakTopologicalOrdering();
    QNodeID start = 0;
    QNodeID end = 0;
    for (QNodeID u = 0; u < qg.numberOfNodes(); ++u) {
      if (ordering[u] == 0) {
        end = u;
      } else if (ordering[u] > ordering[start]) {
        start = u;
      }
    }

    ASSERT_THAT(ordering[end], Eq(0));
    ASSERT_THAT(ordering[start], Gt(0));
    ASSERT_THAT_MOVES_ALONG_PATH_ARE_LEGAL(hg, qg, start, end);
  }

  void ASSERT_THAT_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_RANDOMIZED(const std::string& filename) {
    hg = loadHypergraph(filename);
    applySingletonPartition();
    auto qg = createQuotientGraph();
    qg.reduceToOneRoot();

    QNodeID start = Randomize::instance().getRandomInt(0, qg.numberOfNodes() - 1);
    QNodeID end;
    do {
      end = Randomize::instance().getRandomInt(0, qg.numberOfNodes() - 1);
    } while (end == start);

    ASSERT_THAT_MOVES_ALONG_PATH_ARE_LEGAL(hg, qg, start, end);
  }

  void ASSERT_THAT_MULTIPLE_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_RANDOMIZED(const std::string& filename, PartitionID k) {
    hg = loadHypergraph(filename);
    partitionUsingTopologicalOrdering(k);
    auto qg = createQuotientGraph();
    qg.addMissingEdges();
    qg.reduceToOneRoot();

    QNodeID start;
    QNodeID end;
    std::tie(start, end) = getRandomStartEndInQG(qg);
    ASSERT_THAT_MULTIPLE_MOVES_ALONG_PATH_ARE_LEGAL(hg, qg, start, end);
  }

  void ASSERT_THAT_MULTIPLE_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_EXTENSIVE(const std::string& filename, PartitionID k) {
    hg = loadHypergraph(filename);

    for (QNodeID u = 0; u < k; ++u) {
      for (QNodeID v = 0; v < k; ++v) {
        LOG << V(u) << V(v);
        if (u == v) {
          continue;
        }
        partitionUsingTopologicalOrdering(k);
        auto qg = createQuotientGraph();
        qg.addMissingEdges();
        qg.reduceToOneRoot();
        ASSERT_THAT_QMG_IS_SCC(qg.computeMoveGraph());
        ASSERT_THAT_MULTIPLE_MOVES_ALONG_PATH_ARE_LEGAL(hg, qg, u, v);
      }
    }
  }

  void ASSERT_THAT_MOVES_ALONG_MOVE_GRAPH_EDGES_ARE_LEGAL(const std::string& filename) {
    hg = loadHypergraph(filename);
    applySingletonPartition();
    auto qg = createQuotientGraph();
    qg.reduceToOneRoot();
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

  void ASSERT_THAT_QMG_IS_SCC_SINGLETON_KTS(const std::string& filename, const std::vector<PartitionID> ks) {
    hg = loadHypergraph(filename);
    {
      applySingletonPartition();
      auto qg = createQuotientGraph();
      qg.reduceToOneRoot();
      const auto qmg = qg.computeMoveGraph();
      ASSERT_THAT_QMG_IS_SCC(qmg);
    }
    for (const PartitionID& k : ks) {
      hg.resetPartitioning();
      partitionUsingTopologicalOrdering(k);
      auto qg = createQuotientGraph();
      qg.reduceToOneRoot();
      const auto qmg = qg.computeMoveGraph();
      ASSERT_THAT_QMG_IS_SCC(qmg);
    }
  }

 private:
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

  static void ASSERT_THAT_MOVES_ALONG_PATH_ARE_LEGAL(Hypergraph& hg, QuotientGraph<DFSCycleDetector>& qg,
                                                     QNodeID start, QNodeID end) {
    const auto qmg = qg.computeMoveGraph();
    std::vector<QNodeID> path;
    bool found_path = qmg.findPath(start, end, path);
    ASSERT_THAT(found_path, Eq(true));

    for (std::size_t i = 1; i < path.size(); ++i) {
      QNodeID from = path[i - 1];
      QNodeID to = path[i];
      bool ok = qg.update(from, from, to);
      ASSERT_THAT(ok, Eq(true));
      hg.changeNodePart(from, from, to);
    }
  }

  static void ASSERT_THAT_MULTIPLE_MOVES_ALONG_PATH_ARE_LEGAL(Hypergraph& hg, QuotientGraph<DFSCycleDetector>& qg,
                                                              const QNodeID start, const QNodeID end) {
    const auto qmg = qg.computeMoveGraph();
    const auto ordering = qg.computeWeakTopologicalOrdering();

    std::vector<QNodeID> path;
    bool found_path = qmg.findPath(start, end, path);
    ASSERT_THAT(found_path, Eq(true)) << V(start) << " " << V(end);

    for (std::size_t i = path.size() - 1; i > 0; --i) {
      const QNodeID from = path[i - 1];
      const QNodeID to = path[i];
      movePartitionHalf(hg, qg, from, to, ordering[from] < ordering[to]);
    }
  }

  static void movePartitionHalf(Hypergraph& hg, QuotientGraph<DFSCycleDetector>& qg,
                                const PartitionID from, const PartitionID to, const bool lower) {
    std::size_t count = 0;
    std::size_t limit = hg.partWeight(from) / 2;

    while (count < limit / 2) {
      // find candidate
      HypernodeID candidate = hg.currentNumNodes();
      for (const HypernodeID& hn : hg.nodes()) {
        if (hg.partID(hn) != from) {
          continue;
        }
        bool movable = true;
        if (lower) {
          for (const HyperedgeID& he : hg.incidentTailEdges(hn)) {
            for (const HypernodeID& ht : hg.heads(he)) {
              if (hg.partID(ht) == from) {
                movable = false;
                break;
              }
            }
            if (!movable) {
              break;
            }
          }
        } else {
          for (const HyperedgeID& he : hg.incidentHeadEdges(hn)) {
            for (const HypernodeID& hh : hg.tails(he)) {
              if (hg.partID(hh) == from) {
                movable = false;
                break;
              }
            }
            if (!movable) {
              break;
            }
          }
        }

        if (movable) {
          candidate = hn;
          break;
        }
      }
      ASSERT_THAT(candidate, Lt(hg.currentNumNodes()))
                << "not enough candidates: " << V(from) << " " << V(to) << " " << V(count) << " " << V(limit);
      bool move_ok = qg.update(candidate, from, to);
      ASSERT_THAT(move_ok, Eq(true))
                << V(candidate) << " " << V(from) << " " << V(to) << " " << V(count) << " "
                << V(hg.partID(candidate));
      hg.changeNodePart(candidate, from, to);
      ++count;
    }
  }

  static std::pair<QNodeID, QNodeID> getRandomStartEndInQG(const QuotientGraph<DFSCycleDetector>& qg) {
    QNodeID start = Randomize::instance().getRandomInt(0, qg.numberOfNodes() - 1);
    QNodeID end;
    do {
      end = Randomize::instance().getRandomInt(0, qg.numberOfNodes() - 1);
    } while (end == start);
    return {start, end};
  }

 protected:
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

TEST_F(QuotientGraphTest, C17_MovesAlongMoveGraphAreLegal) {
  ASSERT_THAT_MOVES_ALONG_MOVE_GRAPH_EDGES_ARE_LEGAL("test_instances/c17.hgr");
}

TEST_F(QuotientGraphTest, C17_MultipleMovesAlongPathInMoveGraphAreLegal) {
  constexpr std::size_t N = 10;
  for (std::size_t i = 0; i < N; ++i) {
    ASSERT_THAT_MULTIPLE_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_RANDOMIZED("test_instances/c17.hgr", 4);
  }
}

TEST_F(QuotientGraphTest, C17_MovesAlongPathInMoveGraphAreLegal) {
  ASSERT_THAT_MOVES_ALONG_MOVE_PATH_ARE_LEGAL("test_instances/c17.hgr");
}

TEST_F(QuotientGraphTest, C17_MoveGraphIsSCC) {
  ASSERT_THAT_QMG_IS_SCC_SINGLETON_KTS("test_instances/c17.hgr", {2, 4, 8});
}

TEST_F(QuotientGraphTest, C3540_MovesAlongMoveGraphAreLegal) {
  ASSERT_THAT_MOVES_ALONG_MOVE_GRAPH_EDGES_ARE_LEGAL("test_instances/c3540.hgr");
}

TEST_F(QuotientGraphTest, C3540_MovesAlongPathInMoveGraphAreLegal) {
  ASSERT_THAT_MOVES_ALONG_MOVE_PATH_ARE_LEGAL("test_instances/c3540.hgr");
}

TEST_F(QuotientGraphTest, C3540_MoveGraphIsSCC) {
  ASSERT_THAT_QMG_IS_SCC_SINGLETON_KTS("test_instances/c3540.hgr", {2, 4, 8, 16, 32, 64, 128, 256});
}

TEST_F(QuotientGraphTest, C3540_MovesAlongPathInMoveGraphAreLegalRandomized) {
  constexpr std::size_t N = 10;
  for (std::size_t i = 0; i < N; ++i) {
    ASSERT_THAT_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_RANDOMIZED("test_instances/c3540.hgr");
  }
}

TEST_F(QuotientGraphTest, C880_MultipleMovesAlongPathInMoveGraphAreLegal_Extensive) {
  ASSERT_THAT_MULTIPLE_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_EXTENSIVE("test_instances/c880.hgr", 2);
  ASSERT_THAT_MULTIPLE_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_EXTENSIVE("test_instances/c880.hgr", 4);
  ASSERT_THAT_MULTIPLE_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_EXTENSIVE("test_instances/c880.hgr", 8);
  ASSERT_THAT_MULTIPLE_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_EXTENSIVE("test_instances/c880.hgr", 16);
  ASSERT_THAT_MULTIPLE_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_EXTENSIVE("test_instances/c880.hgr", 32);
  ASSERT_THAT_MULTIPLE_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_EXTENSIVE("test_instances/c880.hgr", 64);
}

TEST_F(QuotientGraphTest, C880_MultipleMovesAlongPathInMoveGraphAreLegal) {
  constexpr std::size_t N = 10;
  for (std::size_t i = 0; i < N; ++i) {
    ASSERT_THAT_MULTIPLE_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_RANDOMIZED("test_instances/c880.hgr", 2);
    ASSERT_THAT_MULTIPLE_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_RANDOMIZED("test_instances/c880.hgr", 4);
    ASSERT_THAT_MULTIPLE_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_RANDOMIZED("test_instances/c880.hgr", 8);
    ASSERT_THAT_MULTIPLE_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_RANDOMIZED("test_instances/c880.hgr", 16);
    ASSERT_THAT_MULTIPLE_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_RANDOMIZED("test_instances/c880.hgr", 32);
    ASSERT_THAT_MULTIPLE_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_RANDOMIZED("test_instances/c880.hgr", 64);
    ASSERT_THAT_MULTIPLE_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_RANDOMIZED("test_instances/c880.hgr", 128);
    ASSERT_THAT_MULTIPLE_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_RANDOMIZED("test_instances/c880.hgr", 256);
    ASSERT_THAT_MULTIPLE_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_RANDOMIZED("test_instances/c880.hgr", 512);
  }
}

TEST_F(QuotientGraphTest, C3540_MultipleMovesAlongPathInMoveGraphAreLegal) {
  constexpr std::size_t N = 10;
  for (std::size_t i = 0; i < N; ++i) {
    ASSERT_THAT_MULTIPLE_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_RANDOMIZED("test_instances/c3540.hgr", 2);
    ASSERT_THAT_MULTIPLE_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_RANDOMIZED("test_instances/c3540.hgr", 4);
    ASSERT_THAT_MULTIPLE_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_RANDOMIZED("test_instances/c3540.hgr", 8);
    ASSERT_THAT_MULTIPLE_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_RANDOMIZED("test_instances/c3540.hgr", 16);
    ASSERT_THAT_MULTIPLE_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_RANDOMIZED("test_instances/c3540.hgr", 32);
    ASSERT_THAT_MULTIPLE_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_RANDOMIZED("test_instances/c3540.hgr", 64);
    ASSERT_THAT_MULTIPLE_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_RANDOMIZED("test_instances/c3540.hgr", 128);
    ASSERT_THAT_MULTIPLE_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_RANDOMIZED("test_instances/c3540.hgr", 256);
    ASSERT_THAT_MULTIPLE_MOVES_ALONG_MOVE_PATH_ARE_LEGAL_RANDOMIZED("test_instances/c3540.hgr", 512);
  }
}
} // namespace dag
} // namespace kahypar