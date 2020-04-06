#include "gmock/gmock.h"

#include "kahypar/partition/dhgp/topord.h"
#include "kahypar/definitions.h"

#include "test_instances.h"

using ::testing::Test;
using ::testing::Not;
using ::testing::Contains;
using ::testing::SizeIs;
using ::testing::Lt;

namespace kahypar::dhgp {
class ADirectedAcyclicHypergraph : public Test {
  public:
  ADirectedAcyclicHypergraph() : hg(test_instances::c17()) {}
  Hypergraph hg;
};

TEST_F(ADirectedAcyclicHypergraph, TopOrdIsTopological) {
  const auto topord = calculateTopologicalOrdering(hg);
  ASSERT_THAT(topord, SizeIs(hg.initialNumNodes()));

  std::set<HypernodeID> already_seen;
  for (const HypernodeID& u : topord) {
    for (const HyperedgeID& he : hg.incidentTailEdges(u)) {
      for (const HypernodeID& v : hg.heads(he)) {
        ASSERT_THAT(already_seen, Not(Contains(v)));
      }
    }
    already_seen.insert(u);
  }

  ASSERT_THAT(already_seen, SizeIs(hg.initialNumNodes()));
}

TEST_F(ADirectedAcyclicHypergraph, InvertedTopOrdIsTopological) {
  const auto inverted_topord = calculateInvertedTopologicalOrdering(hg);
  ASSERT_THAT(inverted_topord, SizeIs(hg.initialNumNodes()));

  for (const HypernodeID& u : inverted_topord) {
    for (const HyperedgeID& he : hg.incidentTailEdges(u)) {
      for (const HypernodeID& v : hg.heads(he)) {
        ASSERT_THAT(inverted_topord[u], Lt(inverted_topord[v]));
      }
    }
  }
}
}