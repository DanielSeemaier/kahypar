#include "gmock/gmock.h"

#include "kahypar/definitions.h"
#include "kahypar/partition/context.h"
#include "kahypar/dag/topological_ordering.h"
#include "kahypar/io/partitioning_output.h"

#include "dag.h"

using ::testing::Eq;
using ::testing::Test;
using ::testing::UnorderedElementsAre;
using ::testing::IsEmpty;

namespace kahypar {
namespace dag {
class DirectedContractionTest : public Test {
 protected:
  void SetUp() override {
    hg = loadHypergraph("test_instances/c17.hgr");
    placeAllHypernodesInPartition(hg, 0);
    assertGraphRestored();
  }

  void partitionUsingTopologicalOrdering(const PartitionID k) {
    Randomize::instance().setSeed(0);
    auto ordering = calculateTopologicalOrdering(hg);
    PartitionID part = 0;
    HypernodeID nodes_per_part = hg.currentNumPins() / k;
    HypernodeID nodes_in_cur_part = 0;

    hg.resetPartitioning();
    hg.changeK(k);
    for (const HypernodeID& hn : ordering) {
      hg.setNodePart(hn, part);
      if (nodes_in_cur_part == nodes_per_part) {
        ++part;
        nodes_in_cur_part = 0;
      }
    }
  }

  void assertGraphRestored() const {
    // edges contain the correct heads and tails
    ASSERT_THAT(toVec(hg.heads(0)), UnorderedElementsAre(0));
    ASSERT_THAT(toVec(hg.tails(0)), UnorderedElementsAre(2, 7));
    ASSERT_THAT(toVec(hg.heads(1)), UnorderedElementsAre(1));
    ASSERT_THAT(toVec(hg.tails(1)), UnorderedElementsAre(8, 2));
    ASSERT_THAT(toVec(hg.tails(2)), UnorderedElementsAre(10, 4));
    ASSERT_THAT(toVec(hg.heads(2)), UnorderedElementsAre(2));
    ASSERT_THAT(toVec(hg.heads(3)), UnorderedElementsAre(3));
    ASSERT_THAT(toVec(hg.tails(3)), UnorderedElementsAre(5, 1));
    ASSERT_THAT(toVec(hg.heads(4)), UnorderedElementsAre(5));
    ASSERT_THAT(toVec(hg.tails(4)), UnorderedElementsAre(6, 10));
    ASSERT_THAT(toVec(hg.heads(5)), UnorderedElementsAre(9));
    ASSERT_THAT(toVec(hg.tails(5)), UnorderedElementsAre(1, 0));

    // nodes contain the correct head and tail edges
    ASSERT_THAT(toVec(hg.incidentHeadEdges(0)), UnorderedElementsAre(0));
    ASSERT_THAT(toVec(hg.incidentTailEdges(0)), UnorderedElementsAre(5));
    ASSERT_THAT(toVec(hg.incidentHeadEdges(1)), UnorderedElementsAre(1));
    ASSERT_THAT(toVec(hg.incidentTailEdges(1)), UnorderedElementsAre(3, 5));
    ASSERT_THAT(toVec(hg.incidentHeadEdges(2)), UnorderedElementsAre(2));
    ASSERT_THAT(toVec(hg.incidentTailEdges(2)), UnorderedElementsAre(0, 1));
    ASSERT_THAT(toVec(hg.incidentHeadEdges(3)), UnorderedElementsAre(3));
    ASSERT_THAT(toVec(hg.incidentTailEdges(3)), IsEmpty());
    ASSERT_THAT(toVec(hg.incidentHeadEdges(4)), IsEmpty());
    ASSERT_THAT(toVec(hg.incidentTailEdges(4)), UnorderedElementsAre(2));
    ASSERT_THAT(toVec(hg.incidentHeadEdges(5)), UnorderedElementsAre(4));
    ASSERT_THAT(toVec(hg.incidentTailEdges(5)), UnorderedElementsAre(3));
    ASSERT_THAT(toVec(hg.incidentHeadEdges(6)), IsEmpty());
    ASSERT_THAT(toVec(hg.incidentTailEdges(6)), UnorderedElementsAre(4));
    ASSERT_THAT(toVec(hg.incidentHeadEdges(7)), IsEmpty());
    ASSERT_THAT(toVec(hg.incidentTailEdges(7)), UnorderedElementsAre(0));
    ASSERT_THAT(toVec(hg.incidentHeadEdges(8)), IsEmpty());
    ASSERT_THAT(toVec(hg.incidentTailEdges(8)), UnorderedElementsAre(1));
    ASSERT_THAT(toVec(hg.incidentHeadEdges(9)), UnorderedElementsAre(5));
    ASSERT_THAT(toVec(hg.incidentTailEdges(9)), IsEmpty());
    ASSERT_THAT(toVec(hg.incidentHeadEdges(10)), IsEmpty());
    ASSERT_THAT(toVec(hg.incidentTailEdges(10)), UnorderedElementsAre(2, 4));
  }

  Hypergraph hg;
};

TEST_F(DirectedContractionTest, Simple_TailToTail_Case1_Contraction) {
  const auto memento_2_7 = hg.contract(2, 7);
  ASSERT_THAT(toVec(hg.heads(0)), UnorderedElementsAre(0));
  ASSERT_THAT(toVec(hg.tails(0)), UnorderedElementsAre(2));
  ASSERT_THAT(toVec(hg.tails(1)), UnorderedElementsAre(2, 8));
  ASSERT_THAT(toVec(hg.incidentHeadEdges(2)), UnorderedElementsAre(2));
  ASSERT_THAT(toVec(hg.incidentTailEdges(2)), UnorderedElementsAre(0, 1));
  hg.uncontract(memento_2_7);
  assertGraphRestored();

  const auto memento_7_2 = hg.contract(7, 2);
  ASSERT_THAT(toVec(hg.heads(0)), UnorderedElementsAre(0));
  ASSERT_THAT(toVec(hg.tails(0)), UnorderedElementsAre(7));
  ASSERT_THAT(toVec(hg.tails(1)), UnorderedElementsAre(7, 8));
  ASSERT_THAT(toVec(hg.heads(2)), UnorderedElementsAre(7));
  hg.uncontract(memento_7_2);
  assertGraphRestored();
}

TEST_F(DirectedContractionTest, Simple_TailToTail_Case2_Contraction) {
  const auto memento_1_10 = hg.contract(1, 10);
  ASSERT_THAT(toVec(hg.heads(1)), UnorderedElementsAre(1));
  ASSERT_THAT(toVec(hg.tails(1)), UnorderedElementsAre(8, 2));
  ASSERT_THAT(toVec(hg.heads(2)), UnorderedElementsAre(2));
  ASSERT_THAT(toVec(hg.tails(2)), UnorderedElementsAre(1, 4));
  ASSERT_THAT(toVec(hg.heads(3)), UnorderedElementsAre(3));
  ASSERT_THAT(toVec(hg.tails(3)), UnorderedElementsAre(5, 1));
  ASSERT_THAT(toVec(hg.heads(4)), UnorderedElementsAre(5));
  ASSERT_THAT(toVec(hg.tails(4)), UnorderedElementsAre(6, 1));
  ASSERT_THAT(toVec(hg.heads(5)), UnorderedElementsAre(9));
  ASSERT_THAT(toVec(hg.tails(5)), UnorderedElementsAre(1, 0));
  ASSERT_THAT(toVec(hg.incidentHeadEdges(1)), UnorderedElementsAre(1));
  ASSERT_THAT(toVec(hg.incidentTailEdges(1)), UnorderedElementsAre(2, 3, 4, 5));
  hg.uncontract(memento_1_10);
  assertGraphRestored();
}

TEST_F(DirectedContractionTest, Simple_TailToHead_Case1_Contraction) {
  const auto memento_2_0 = hg.contract(2, 0);
  ASSERT_THAT(toVec(hg.heads(0)), UnorderedElementsAre(2));
  ASSERT_THAT(toVec(hg.tails(0)), UnorderedElementsAre(7));
  ASSERT_THAT(toVec(hg.heads(5)), UnorderedElementsAre(9));
  ASSERT_THAT(toVec(hg.tails(5)), UnorderedElementsAre(1, 2));
  ASSERT_THAT(toVec(hg.incidentHeadEdges(2)), UnorderedElementsAre(0, 2));
  ASSERT_THAT(toVec(hg.incidentTailEdges(2)), UnorderedElementsAre(1, 5));
  hg.uncontract(memento_2_0);
  assertGraphRestored();
}

TEST_F(DirectedContractionTest, Simple_HeadToTail_Case1_Contraction) {
  const auto memento_0_2 = hg.contract(0, 2);
  ASSERT_THAT(toVec(hg.heads(0)), UnorderedElementsAre(0));
  ASSERT_THAT(toVec(hg.tails(0)), UnorderedElementsAre(7));
  ASSERT_THAT(toVec(hg.heads(1)), UnorderedElementsAre(1));
  ASSERT_THAT(toVec(hg.tails(1)), UnorderedElementsAre(8, 0));
  ASSERT_THAT(toVec(hg.heads(2)), UnorderedElementsAre(0));
  ASSERT_THAT(toVec(hg.tails(2)), UnorderedElementsAre(10, 4));
  ASSERT_THAT(toVec(hg.heads(5)), UnorderedElementsAre(9));
  ASSERT_THAT(toVec(hg.tails(5)), UnorderedElementsAre(1, 0));
  ASSERT_THAT(toVec(hg.incidentHeadEdges(0)), UnorderedElementsAre(0, 2));
  ASSERT_THAT(toVec(hg.incidentTailEdges(0)), UnorderedElementsAre(1, 5));
  hg.uncontract(memento_0_2);
  assertGraphRestored();
}

TEST_F(DirectedContractionTest, ContractEdge) {
  // (0, (2, 7))
  const auto memento_l2_7r = hg.contract(2, 7);
  const auto memento_l0_l2_7rr = hg.contract(0, 2);
  ASSERT_THAT(toVec(hg.heads(0)), UnorderedElementsAre(0));
  ASSERT_THAT(toVec(hg.tails(0)), IsEmpty());
  ASSERT_THAT(toVec(hg.heads(1)), UnorderedElementsAre(1));
  ASSERT_THAT(toVec(hg.tails(1)), UnorderedElementsAre(8, 0));
  ASSERT_THAT(toVec(hg.heads(2)), UnorderedElementsAre(0));
  ASSERT_THAT(toVec(hg.tails(2)), UnorderedElementsAre(10, 4));
  ASSERT_THAT(toVec(hg.heads(5)), UnorderedElementsAre(9));
  ASSERT_THAT(toVec(hg.tails(5)), UnorderedElementsAre(1, 0));
  ASSERT_THAT(toVec(hg.incidentHeadEdges(0)), UnorderedElementsAre(0, 2));
  ASSERT_THAT(toVec(hg.incidentTailEdges(0)), UnorderedElementsAre(1, 5));
  hg.uncontract(memento_l0_l2_7rr);
  hg.uncontract(memento_l2_7r);
  assertGraphRestored();

  // ((0, 2), 7)
  const auto memento_l0_2r = hg.contract(0, 2);
  const auto memento_ll0_2r_7r = hg.contract(0, 7);
  ASSERT_THAT(toVec(hg.heads(0)), UnorderedElementsAre(0));
  ASSERT_THAT(toVec(hg.tails(0)), IsEmpty());
  ASSERT_THAT(toVec(hg.heads(1)), UnorderedElementsAre(1));
  ASSERT_THAT(toVec(hg.tails(1)), UnorderedElementsAre(8, 0));
  ASSERT_THAT(toVec(hg.heads(2)), UnorderedElementsAre(0));
  ASSERT_THAT(toVec(hg.tails(2)), UnorderedElementsAre(10, 4));
  ASSERT_THAT(toVec(hg.heads(5)), UnorderedElementsAre(9));
  ASSERT_THAT(toVec(hg.tails(5)), UnorderedElementsAre(1, 0));
  ASSERT_THAT(toVec(hg.incidentHeadEdges(0)), UnorderedElementsAre(0, 2));
  ASSERT_THAT(toVec(hg.incidentTailEdges(0)), UnorderedElementsAre(1, 5));
  hg.uncontract(memento_ll0_2r_7r);
  hg.uncontract(memento_l0_2r);
  assertGraphRestored();
}

TEST_F(DirectedContractionTest, ContractGraph) {
  std::vector<Hypergraph::Memento> mementos{
    hg.contract(0, 1),
    hg.contract(0, 2),
    hg.contract(0, 3),
    hg.contract(0, 4),
    hg.contract(0, 5),
    hg.contract(0, 6),
    hg.contract(0, 7),
    hg.contract(0, 8),
    hg.contract(0, 9),
    hg.contract(0, 10),
  };
  ASSERT_THAT(toVec(hg.heads(0)), UnorderedElementsAre(0));
  ASSERT_THAT(toVec(hg.tails(0)), IsEmpty());
  ASSERT_THAT(toVec(hg.heads(1)), UnorderedElementsAre(0));
  ASSERT_THAT(toVec(hg.tails(1)), IsEmpty());
  ASSERT_THAT(toVec(hg.heads(2)), UnorderedElementsAre(0));
  ASSERT_THAT(toVec(hg.tails(2)), IsEmpty());
  ASSERT_THAT(toVec(hg.heads(3)), UnorderedElementsAre(0));
  ASSERT_THAT(toVec(hg.tails(3)), IsEmpty());
  ASSERT_THAT(toVec(hg.heads(4)), UnorderedElementsAre(0));
  ASSERT_THAT(toVec(hg.tails(4)), IsEmpty());
  ASSERT_THAT(toVec(hg.heads(5)), UnorderedElementsAre(0));
  ASSERT_THAT(toVec(hg.tails(5)), IsEmpty());
  ASSERT_THAT(toVec(hg.incidentHeadEdges(0)), UnorderedElementsAre(0, 1, 2, 3, 4, 5));
  ASSERT_THAT(toVec(hg.incidentTailEdges(0)), IsEmpty());
  for (std::size_t i = mementos.size(); i > 0; --i) {
    hg.uncontract(mementos[i - 1]);
  }
  assertGraphRestored();
}

TEST_F(DirectedContractionTest, RemoveAndRestoreEdge) {
  for (std::size_t i = 0; i < 2; ++i) {
    hg.removeEdge(0);
    ASSERT_THAT(hg.currentNumEdges(), Eq(5));
    ASSERT_THAT(toVec(hg.incidentHeadEdges(0)), IsEmpty());
    ASSERT_THAT(toVec(hg.incidentTailEdges(2)), UnorderedElementsAre(1));
    ASSERT_THAT(toVec(hg.incidentTailEdges(7)), IsEmpty());
    hg.restoreEdge(0);
    assertGraphRestored();
  }
  for (std::size_t i = 0; i < 2; ++i) {
    hg.removeEdge(0);
    hg.removeEdge(1);
    hg.removeEdge(2);
    ASSERT_THAT(hg.currentNumEdges(), Eq(3));
    ASSERT_THAT(toVec(hg.incidentHeadEdges(0)), IsEmpty());
    ASSERT_THAT(toVec(hg.incidentTailEdges(0)), UnorderedElementsAre(5));
    ASSERT_THAT(toVec(hg.incidentHeadEdges(1)), IsEmpty());
    ASSERT_THAT(toVec(hg.incidentTailEdges(1)), UnorderedElementsAre(3, 5));
    ASSERT_THAT(toVec(hg.incidentHeadEdges(2)), IsEmpty());
    ASSERT_THAT(toVec(hg.incidentTailEdges(2)), IsEmpty());
    ASSERT_THAT(toVec(hg.incidentHeadEdges(3)), UnorderedElementsAre(3));
    ASSERT_THAT(toVec(hg.incidentTailEdges(3)), IsEmpty());
    ASSERT_THAT(toVec(hg.incidentHeadEdges(4)), IsEmpty());
    ASSERT_THAT(toVec(hg.incidentTailEdges(4)), IsEmpty());
    ASSERT_THAT(toVec(hg.incidentHeadEdges(5)), UnorderedElementsAre(4));
    ASSERT_THAT(toVec(hg.incidentTailEdges(5)), UnorderedElementsAre(3));
    ASSERT_THAT(toVec(hg.incidentHeadEdges(6)), IsEmpty());
    ASSERT_THAT(toVec(hg.incidentTailEdges(6)), UnorderedElementsAre(4));
    ASSERT_THAT(toVec(hg.incidentHeadEdges(7)), IsEmpty());
    ASSERT_THAT(toVec(hg.incidentTailEdges(7)), IsEmpty());
    ASSERT_THAT(toVec(hg.incidentHeadEdges(8)), IsEmpty());
    ASSERT_THAT(toVec(hg.incidentTailEdges(8)), IsEmpty());
    ASSERT_THAT(toVec(hg.incidentHeadEdges(9)), UnorderedElementsAre(5));
    ASSERT_THAT(toVec(hg.incidentTailEdges(9)), IsEmpty());
    ASSERT_THAT(toVec(hg.incidentHeadEdges(10)), IsEmpty());
    ASSERT_THAT(toVec(hg.incidentTailEdges(10)), UnorderedElementsAre(4));
    hg.restoreEdge(2);
    hg.restoreEdge(1);
    hg.restoreEdge(0);
    assertGraphRestored();
  }
}

TEST_F(DirectedContractionTest, RemoveRestoreSingletonHypernodeInCoarserGraph) {
  auto memento_1 = hg.contract(2, 7);
  auto memento_2 = hg.contract(0, 2);
  hg.removeEdge(0);
  ASSERT_THAT(hg.currentNumEdges(), Eq(5));
  ASSERT_THAT(toVec(hg.incidentEdges(0)), UnorderedElementsAre(1, 2, 5));
  ASSERT_THAT(toVec(hg.incidentHeadEdges(0)), UnorderedElementsAre(2));
  ASSERT_THAT(toVec(hg.incidentTailEdges(0)), UnorderedElementsAre(1, 5));
  ASSERT_THAT(toVec(hg.heads(2)), UnorderedElementsAre(0));
  ASSERT_THAT(toVec(hg.tails(2)), UnorderedElementsAre(10, 4));
  hg.restoreEdge(0);
  ASSERT_THAT(hg.currentNumEdges(), Eq(6));
  ASSERT_THAT(toVec(hg.incidentEdges(0)), UnorderedElementsAre(0, 1, 2, 5));
  ASSERT_THAT(toVec(hg.incidentHeadEdges(0)), UnorderedElementsAre(0, 2));
  ASSERT_THAT(toVec(hg.incidentTailEdges(0)), UnorderedElementsAre(1, 5));
  ASSERT_THAT(toVec(hg.pins(2)), UnorderedElementsAre(0, 4, 10));
  ASSERT_THAT(toVec(hg.heads(2)), UnorderedElementsAre(0));
  ASSERT_THAT(toVec(hg.tails(2)), UnorderedElementsAre(10, 4));
  hg.uncontract(memento_2);
  hg.uncontract(memento_1);
}

// Code that crashed at some point in time

TEST_F(DirectedContractionTest, Regression_C3540_UncoarseningCrash) {
  hg = loadHypergraph("test_instances/c3540.hgr");
  partitionUsingTopologicalOrdering(8);

  auto memento = hg.contract(1604, 1566);
  hg.removeEdge(1565);
  hg.initializeNumCutHyperedges();
  hg.restoreEdge(1565);
  hg.uncontract(memento); // should not crash
}
} // namespace dag
} // namespace kahypar
