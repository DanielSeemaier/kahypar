#include "gmock/gmock.h"

#include "dag.h"

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
} // namespace dag
} // namespace kahypar
