#include "gmock/gmock.h"

#include "kahypar/definitions.h"

#include "../datastructure/hypergraph_test_fixtures.h"
#include "matcher.h"
#include "test_instances.h"

using ::testing::Eq;
using ::testing::IsEmpty;
using ::testing::Test;
using ::testing::UnorderedElementsAre;

namespace kahypar::ds {
using namespace kahypar::dhgp;
using testing::IsSameDirectedHypergraph;
using testing::utility::iter2vec;

class ADirectedHypergraph : public Test {
 public:
  ADirectedHypergraph() : hypergraph(test_instances::c17()) {}
  Hypergraph hypergraph;
};

TEST_F(ADirectedHypergraph, DirectedHypergraphIsDirected) {
  ASSERT_TRUE(hypergraph.isDirected());
}

TEST_F(AHypergraph, UndirectedHypergraphIsUndirected) {
  ASSERT_FALSE(hypergraph.isDirected());
}

TEST_F(ADirectedHypergraph, HasCorrectNumberOfHeadsAndTails) {
  for (const HyperedgeID he : hypergraph.edges()) {
    ASSERT_THAT(hypergraph.edgeNumHeadPins(he), Eq(1));
    ASSERT_THAT(hypergraph.edgeNumTailPins(he), Eq(2));
  }
  for (const auto& hn : {0, 1, 2, 3, 5, 9}) {
    ASSERT_THAT(hypergraph.nodeNumHeadEdges(hn), Eq(1));
  }
  for (const auto& hn : {1, 2, 10}) {
    ASSERT_THAT(hypergraph.nodeNumTailEdges(hn), Eq(2));
  }
  for (const auto& hn : {0, 4, 5, 6, 7, 8}) {
    ASSERT_THAT(hypergraph.nodeNumTailEdges(hn), Eq(1));
  }
}

TEST_F(ADirectedHypergraph, StructureIsCorrect) {
  // edges contain the correct heads and tails
  ASSERT_THAT(iter2vec(hypergraph.headPins(0)), UnorderedElementsAre(0));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(0)), UnorderedElementsAre(2, 7));
  ASSERT_THAT(iter2vec(hypergraph.headPins(1)), UnorderedElementsAre(1));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(1)), UnorderedElementsAre(8, 2));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(2)), UnorderedElementsAre(10, 4));
  ASSERT_THAT(iter2vec(hypergraph.headPins(2)), UnorderedElementsAre(2));
  ASSERT_THAT(iter2vec(hypergraph.headPins(3)), UnorderedElementsAre(3));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(3)), UnorderedElementsAre(5, 1));
  ASSERT_THAT(iter2vec(hypergraph.headPins(4)), UnorderedElementsAre(5));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(4)), UnorderedElementsAre(6, 10));
  ASSERT_THAT(iter2vec(hypergraph.headPins(5)), UnorderedElementsAre(9));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(5)), UnorderedElementsAre(1, 0));

  // nodes contain the correct head and tail edges
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(0)), UnorderedElementsAre(0));
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(0)), UnorderedElementsAre(5));
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(1)), UnorderedElementsAre(1));
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(1)), UnorderedElementsAre(3, 5));
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(2)), UnorderedElementsAre(2));
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(2)), UnorderedElementsAre(0, 1));
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(3)), UnorderedElementsAre(3));
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(3)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(4)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(4)), UnorderedElementsAre(2));
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(5)), UnorderedElementsAre(4));
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(5)), UnorderedElementsAre(3));
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(6)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(6)), UnorderedElementsAre(4));
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(7)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(7)), UnorderedElementsAre(0));
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(8)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(8)), UnorderedElementsAre(1));
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(9)), UnorderedElementsAre(5));
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(9)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(10)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(10)), UnorderedElementsAre(2, 4));
}

TEST_F(ADirectedHypergraph, TailTail_Case1_Contraction) {
  hypergraph.contract(2, 7);
  ASSERT_THAT(iter2vec(hypergraph.headPins(0)), UnorderedElementsAre(0));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(0)), UnorderedElementsAre(2));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(1)), UnorderedElementsAre(2, 8));
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(2)), UnorderedElementsAre(2));
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(2)), UnorderedElementsAre(0, 1));
}

TEST_F(ADirectedHypergraph, TailToTail_Case1_Uncontraction) {
  hypergraph.uncontract(hypergraph.contract(2, 7));
  ASSERT_THAT(hypergraph, IsSameDirectedHypergraph(test_instances::c17()));
}

TEST_F(ADirectedHypergraph, TailTail_Case2_Contraction) {
  hypergraph.contract(1, 10);
  ASSERT_THAT(iter2vec(hypergraph.headPins(1)), UnorderedElementsAre(1));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(1)), UnorderedElementsAre(8, 2));
  ASSERT_THAT(iter2vec(hypergraph.headPins(2)), UnorderedElementsAre(2));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(2)), UnorderedElementsAre(1, 4));
  ASSERT_THAT(iter2vec(hypergraph.headPins(3)), UnorderedElementsAre(3));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(3)), UnorderedElementsAre(5, 1));
  ASSERT_THAT(iter2vec(hypergraph.headPins(4)), UnorderedElementsAre(5));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(4)), UnorderedElementsAre(6, 1));
  ASSERT_THAT(iter2vec(hypergraph.headPins(5)), UnorderedElementsAre(9));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(5)), UnorderedElementsAre(1, 0));
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(1)), UnorderedElementsAre(1));
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(1)), UnorderedElementsAre(2, 3, 4, 5));
}

TEST_F(ADirectedHypergraph, TailTail_Case2_Uncontraction) {
  hypergraph.uncontract(hypergraph.contract(1, 10));
  ASSERT_THAT(hypergraph, IsSameDirectedHypergraph(test_instances::c17()));
}

TEST_F(ADirectedHypergraph, TailHead_Case1_Contraction) {
  hypergraph.contract(2, 0);
  ASSERT_THAT(iter2vec(hypergraph.headPins(0)), UnorderedElementsAre(2));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(0)), UnorderedElementsAre(7));
  ASSERT_THAT(iter2vec(hypergraph.headPins(5)), UnorderedElementsAre(9));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(5)), UnorderedElementsAre(1, 2));
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(2)), UnorderedElementsAre(0, 2));
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(2)), UnorderedElementsAre(1, 5));
}

TEST_F(ADirectedHypergraph, TailHead_Case1_Uncontraction) {
  hypergraph.uncontract(hypergraph.contract(2, 0));
  ASSERT_THAT(hypergraph, IsSameDirectedHypergraph(test_instances::c17()));
}

TEST_F(ADirectedHypergraph, HeadTail_Case1_Contraction) {
  hypergraph.contract(0, 2);
  ASSERT_THAT(iter2vec(hypergraph.headPins(0)), UnorderedElementsAre(0));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(0)), UnorderedElementsAre(7));
  ASSERT_THAT(iter2vec(hypergraph.headPins(1)), UnorderedElementsAre(1));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(1)), UnorderedElementsAre(8, 0));
  ASSERT_THAT(iter2vec(hypergraph.headPins(2)), UnorderedElementsAre(0));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(2)), UnorderedElementsAre(10, 4));
  ASSERT_THAT(iter2vec(hypergraph.headPins(5)), UnorderedElementsAre(9));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(5)), UnorderedElementsAre(1, 0));
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(0)), UnorderedElementsAre(0, 2));
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(0)), UnorderedElementsAre(1, 5));
}

TEST_F(ADirectedHypergraph, HeadTail_Case1_Uncontraction) {
  hypergraph.uncontract(hypergraph.contract(0, 2));
  ASSERT_THAT(hypergraph, IsSameDirectedHypergraph(test_instances::c17()));
}

TEST_F(ADirectedHypergraph, ContractAndUncontractTwoPairs) {
  // (0, (2, 7))
  const auto memento_l2_7r = hypergraph.contract(2, 7);
  const auto memento_l0_l2_7rr = hypergraph.contract(0, 2);
  ASSERT_THAT(iter2vec(hypergraph.headPins(0)), UnorderedElementsAre(0));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(0)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.headPins(1)), UnorderedElementsAre(1));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(1)), UnorderedElementsAre(8, 0));
  ASSERT_THAT(iter2vec(hypergraph.headPins(2)), UnorderedElementsAre(0));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(2)), UnorderedElementsAre(10, 4));
  ASSERT_THAT(iter2vec(hypergraph.headPins(5)), UnorderedElementsAre(9));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(5)), UnorderedElementsAre(1, 0));
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(0)), UnorderedElementsAre(0, 2));
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(0)), UnorderedElementsAre(1, 5));
  hypergraph.uncontract(memento_l0_l2_7rr);
  hypergraph.uncontract(memento_l2_7r);
  ASSERT_THAT(hypergraph, IsSameDirectedHypergraph(test_instances::c17()));

  // ((0, 2), 7)
  const auto memento_l0_2r = hypergraph.contract(0, 2);
  const auto memento_ll0_2r_7r = hypergraph.contract(0, 7);
  ASSERT_THAT(iter2vec(hypergraph.headPins(0)), UnorderedElementsAre(0));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(0)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.headPins(1)), UnorderedElementsAre(1));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(1)), UnorderedElementsAre(8, 0));
  ASSERT_THAT(iter2vec(hypergraph.headPins(2)), UnorderedElementsAre(0));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(2)), UnorderedElementsAre(10, 4));
  ASSERT_THAT(iter2vec(hypergraph.headPins(5)), UnorderedElementsAre(9));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(5)), UnorderedElementsAre(1, 0));
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(0)), UnorderedElementsAre(0, 2));
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(0)), UnorderedElementsAre(1, 5));
  hypergraph.uncontract(memento_ll0_2r_7r);
  hypergraph.uncontract(memento_l0_2r);
  ASSERT_THAT(hypergraph, IsSameDirectedHypergraph(test_instances::c17()));
}

TEST_F(ADirectedHypergraph, ContractAndUncontractGraph) {
  std::vector<Hypergraph::Memento> mementos{
      hypergraph.contract(0, 1),
      hypergraph.contract(0, 2),
      hypergraph.contract(0, 3),
      hypergraph.contract(0, 4),
      hypergraph.contract(0, 5),
      hypergraph.contract(0, 6),
      hypergraph.contract(0, 7),
      hypergraph.contract(0, 8),
      hypergraph.contract(0, 9),
      hypergraph.contract(0, 10),
  };
  ASSERT_THAT(iter2vec(hypergraph.headPins(0)), UnorderedElementsAre(0));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(0)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.headPins(1)), UnorderedElementsAre(0));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(1)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.headPins(2)), UnorderedElementsAre(0));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(2)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.headPins(3)), UnorderedElementsAre(0));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(3)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.headPins(4)), UnorderedElementsAre(0));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(4)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.headPins(5)), UnorderedElementsAre(0));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(5)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(0)), UnorderedElementsAre(0, 1, 2, 3, 4, 5));
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(0)), IsEmpty());
  for (std::size_t i = mementos.size(); i > 0; --i) {
    hypergraph.uncontract(mementos[i - 1]);
  }
  ASSERT_THAT(hypergraph, IsSameDirectedHypergraph(test_instances::c17()));
}

TEST_F(ADirectedHypergraph, RemoveOneEdge) {
  hypergraph.removeEdge(0);
  ASSERT_THAT(hypergraph.currentNumEdges(), Eq(5));
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(0)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(2)), UnorderedElementsAre(1));
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(7)), IsEmpty());
}

TEST_F(ADirectedHypergraph, RestoreOneEdge) {
  hypergraph.removeEdge(0);
  hypergraph.restoreEdge(0);
  ASSERT_THAT(hypergraph, IsSameDirectedHypergraph(test_instances::c17()));
}

TEST_F(ADirectedHypergraph, RemoveAndRestoreTwoEdges) {
  hypergraph.removeEdge(0);
  hypergraph.removeEdge(1);
  hypergraph.removeEdge(2);
  ASSERT_THAT(hypergraph.currentNumEdges(), Eq(3));
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(0)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(0)), UnorderedElementsAre(5));
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(1)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(1)), UnorderedElementsAre(3, 5));
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(2)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(2)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(3)), UnorderedElementsAre(3));
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(3)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(4)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(4)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(5)), UnorderedElementsAre(4));
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(5)), UnorderedElementsAre(3));
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(6)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(6)), UnorderedElementsAre(4));
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(7)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(7)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(8)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(8)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(9)), UnorderedElementsAre(5));
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(9)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(10)), IsEmpty());
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(10)), UnorderedElementsAre(4));
  hypergraph.restoreEdge(2);
  hypergraph.restoreEdge(1);
  hypergraph.restoreEdge(0);
  ASSERT_THAT(hypergraph, IsSameDirectedHypergraph(test_instances::c17()));
}

TEST_F(ADirectedHypergraph, RemoveAndRestoreOneEdgeInCoarserGraph) {
  auto memento_1 = hypergraph.contract(2, 7);
  auto memento_2 = hypergraph.contract(0, 2);
  hypergraph.removeEdge(0);
  ASSERT_THAT(hypergraph.currentNumEdges(), Eq(5));
  ASSERT_THAT(iter2vec(hypergraph.incidentEdges(0)), UnorderedElementsAre(1, 2, 5));
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(0)), UnorderedElementsAre(2));
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(0)), UnorderedElementsAre(1, 5));
  ASSERT_THAT(iter2vec(hypergraph.headPins(2)), UnorderedElementsAre(0));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(2)), UnorderedElementsAre(10, 4));
  hypergraph.restoreEdge(0);
  ASSERT_THAT(hypergraph.currentNumEdges(), Eq(6));
  ASSERT_THAT(iter2vec(hypergraph.incidentEdges(0)), UnorderedElementsAre(0, 1, 2, 5));
  ASSERT_THAT(iter2vec(hypergraph.incidentHeadEdges(0)), UnorderedElementsAre(0, 2));
  ASSERT_THAT(iter2vec(hypergraph.incidentTailEdges(0)), UnorderedElementsAre(1, 5));
  ASSERT_THAT(iter2vec(hypergraph.pins(2)), UnorderedElementsAre(0, 4, 10));
  ASSERT_THAT(iter2vec(hypergraph.headPins(2)), UnorderedElementsAre(0));
  ASSERT_THAT(iter2vec(hypergraph.tailPins(2)), UnorderedElementsAre(10, 4));
  hypergraph.uncontract(memento_2);
  hypergraph.uncontract(memento_1);
  ASSERT_THAT(hypergraph, IsSameDirectedHypergraph(test_instances::c17()));
}
}// namespace kahypar::ds