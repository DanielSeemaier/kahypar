#include "gmock/gmock.h"

#include "kahypar/definitions.h"

#include "../datastructure/hypergraph_test_fixtures.h"

using ::testing::Test;
using ::testing::Eq;
using ::testing::UnorderedElementsAre;
using ::testing::IsEmpty;

namespace kahypar::ds {
class ADirectedHypergraph : public Test {
   public:
  ADirectedHypergraph() :
    hypergraph(11, 6, HyperedgeIndexVector { 0, 3, 6, 9, 12, 15, 18 },
               HyperedgeVector { 0, 2, 7, 1, 8, 2, 2, 10, 4, 3, 5, 1, 5, 6, 10, 9, 1, 0 },
               true, HyperedgeVector { 1, 1, 1, 1, 1, 1 }) { }
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
    ASSERT_THAT(hypergraph.edgeNumHeads(he), Eq(1));
    ASSERT_THAT(hypergraph.edgeNumTails(he), Eq(2));
  }
  for (const auto& hn : {0, 1, 2, 3, 5, 9}) {
    ASSERT_THAT(hypergraph.nodeNumHeads(hn), Eq(1));
  }
  for (const auto& hn : {1, 2, 10}) {
    ASSERT_THAT(hypergraph.nodeNumTails(hn), Eq(2));
  }
  for (const auto& hn : {0, 4, 5, 6, 7, 8}) {
    ASSERT_THAT(hypergraph.nodeNumTails(hn), Eq(1));
  }
}

TEST_F(ADirectedHypergraph, StructureIsCorrect) {
  auto to_vec = [](const auto&& iter_pair) -> auto {
    std::vector<typename decltype(iter_pair.first)::value_type> values;
    values.insert(values.begin(), iter_pair.first, iter_pair.second);
    return values;
  };

  // edges contain the correct heads and tails
  ASSERT_THAT(to_vec(hypergraph.heads(0)), UnorderedElementsAre(0));
  ASSERT_THAT(to_vec(hypergraph.tails(0)), UnorderedElementsAre(2, 7));
  ASSERT_THAT(to_vec(hypergraph.heads(1)), UnorderedElementsAre(1));
  ASSERT_THAT(to_vec(hypergraph.tails(1)), UnorderedElementsAre(8, 2));
  ASSERT_THAT(to_vec(hypergraph.tails(2)), UnorderedElementsAre(10, 4));
  ASSERT_THAT(to_vec(hypergraph.heads(2)), UnorderedElementsAre(2));
  ASSERT_THAT(to_vec(hypergraph.heads(3)), UnorderedElementsAre(3));
  ASSERT_THAT(to_vec(hypergraph.tails(3)), UnorderedElementsAre(5, 1));
  ASSERT_THAT(to_vec(hypergraph.heads(4)), UnorderedElementsAre(5));
  ASSERT_THAT(to_vec(hypergraph.tails(4)), UnorderedElementsAre(6, 10));
  ASSERT_THAT(to_vec(hypergraph.heads(5)), UnorderedElementsAre(9));
  ASSERT_THAT(to_vec(hypergraph.tails(5)), UnorderedElementsAre(1, 0));

  // nodes contain the correct head and tail edges
  ASSERT_THAT(to_vec(hypergraph.incidentHeadEdges(0)), UnorderedElementsAre(0));
  ASSERT_THAT(to_vec(hypergraph.incidentTailEdges(0)), UnorderedElementsAre(5));
  ASSERT_THAT(to_vec(hypergraph.incidentHeadEdges(1)), UnorderedElementsAre(1));
  ASSERT_THAT(to_vec(hypergraph.incidentTailEdges(1)), UnorderedElementsAre(3, 5));
  ASSERT_THAT(to_vec(hypergraph.incidentHeadEdges(2)), UnorderedElementsAre(2));
  ASSERT_THAT(to_vec(hypergraph.incidentTailEdges(2)), UnorderedElementsAre(0, 1));
  ASSERT_THAT(to_vec(hypergraph.incidentHeadEdges(3)), UnorderedElementsAre(3));
  ASSERT_THAT(to_vec(hypergraph.incidentTailEdges(3)), IsEmpty());
  ASSERT_THAT(to_vec(hypergraph.incidentHeadEdges(4)), IsEmpty());
  ASSERT_THAT(to_vec(hypergraph.incidentTailEdges(4)), UnorderedElementsAre(2));
  ASSERT_THAT(to_vec(hypergraph.incidentHeadEdges(5)), UnorderedElementsAre(4));
  ASSERT_THAT(to_vec(hypergraph.incidentTailEdges(5)), UnorderedElementsAre(3));
  ASSERT_THAT(to_vec(hypergraph.incidentHeadEdges(6)), IsEmpty());
  ASSERT_THAT(to_vec(hypergraph.incidentTailEdges(6)), UnorderedElementsAre(4));
  ASSERT_THAT(to_vec(hypergraph.incidentHeadEdges(7)), IsEmpty());
  ASSERT_THAT(to_vec(hypergraph.incidentTailEdges(7)), UnorderedElementsAre(0));
  ASSERT_THAT(to_vec(hypergraph.incidentHeadEdges(8)), IsEmpty());
  ASSERT_THAT(to_vec(hypergraph.incidentTailEdges(8)), UnorderedElementsAre(1));
  ASSERT_THAT(to_vec(hypergraph.incidentHeadEdges(9)), UnorderedElementsAre(5));
  ASSERT_THAT(to_vec(hypergraph.incidentTailEdges(9)), IsEmpty());
  ASSERT_THAT(to_vec(hypergraph.incidentHeadEdges(10)), IsEmpty());
  ASSERT_THAT(to_vec(hypergraph.incidentTailEdges(10)), UnorderedElementsAre(2, 4));
}
} // namespace kahypar::ds