#include "gmock/gmock.h"

#include "kahypar/definitions.h"
#include "kahypar/partition/context.h"
#include "kahypar/dag/topological_ordering.h"
#include "kahypar/io/partitioning_output.h"
#include "kahypar/partition/refinement/km1_gain_manager.h"

#include "dag.h"

using ::testing::Eq;
using ::testing::Ne;
using ::testing::Test;
using ::testing::TestWithParam;
using ::testing::UnorderedElementsAre;
using ::testing::IsEmpty;
using ::testing::Values;
using ::testing::_;

namespace kahypar {
namespace dag {
class KMinusOneGainManagerBaseTest {
 protected:
  void loadGraph(const std::string& filename, const PartitionID k) {
    hg = loadHypergraph(filename);
    partitionUsingTopologicalOrdering(k);
    context.partition.k = k;
    manager = std::make_unique<KMinusOneGainManager>(hg, context);
    manager->initialize();
  }

  void partitionUsingTopologicalOrdering(const PartitionID k) {
    Randomize::instance().setSeed(0);
    auto ordering = calculateTopologicalOrdering(hg);
    PartitionID part = 0;
    HypernodeID nodes_per_part = hg.initialNumNodes() / k + 1;
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

  struct Move {
    const HypernodeID hn;
    const PartitionID from;
    const PartitionID to;

    Move(const HypernodeID hn, const PartitionID from, const PartitionID to) :
      hn(hn),
      from(from),
      to(to) {}
  };

  std::vector<Move> performRandomMoves(const std::size_t num) {
    return performRandomMoves(num,
                              [](const auto&, const auto&, const auto&) {},
                              [](const auto&, const auto&, const auto&) {},
                              [](const auto&, const auto&, const auto&) {},
                              [](const auto&, const auto&, const auto&) {});
  }

  template<typename InsertCallback, typename UpdateCallback, typename RemoveCallback, typename MoveCallback>
  std::vector<Move> performRandomMoves(const std::size_t num, InsertCallback&& insert_cb, UpdateCallback&& update_cb,
                          RemoveCallback&& remove_cb, MoveCallback&& move_cb) {
    std::vector<Move> moves;
    auto& rand = Randomize::instance();
    for (std::size_t i = 0; i < num; ++i) {
      const HypernodeID hn = rand.getRandomInt(0, hg.initialNumNodes() - 1);
      const PartitionID from_part = hg.partID(hn);
      PartitionID to_part;
      do {
        to_part = rand.getRandomInt(0, context.partition.k - 1);
      } while (from_part == to_part);

      hg.changeNodePart(hn, from_part, to_part);
      manager->updateAfterMovement(hn, from_part, to_part, insert_cb, update_cb, remove_cb);
      move_cb(hn, from_part, to_part);
      moves.emplace_back(hn, from_part, to_part);
    }

    return moves;
  }

  Gain gainInducedByHypergraph(const HypernodeID hn, const PartitionID target_part) const {
    Gain gain = 0;
    for (const HyperedgeID& he : hg.incidentEdges(hn)) {
      // our contraction test does not remove single pin HEs, so this assertion might fail
      //ASSERT(hg.edgeSize(he) > 1, V(he));
      gain += gainInducedByHyperedge(hn, he, target_part);
    }
    return gain;
  }

  Gain gainInducedByHyperedge(const HypernodeID hn, const HyperedgeID he, const PartitionID target_part) const {
    const HypernodeID pins_in_source_part = hg.pinCountInPart(he, hg.partID(hn));
    const HypernodeID pins_in_target_part = hg.pinCountInPart(he, target_part);
    const HyperedgeWeight he_weight = hg.edgeWeight(he);
    Gain gain = pins_in_source_part == 1 ? he_weight : 0;
    gain -= pins_in_target_part == 0 ? he_weight : 0;
    return gain;
  }

  bool hypernodeIsConnectedToPart(const HypernodeID hn, const PartitionID part) const {
    for (const HyperedgeID& he : hg.incidentEdges(hn)) {
      if (hg.pinCountInPart(he, part) > 0) {
        return true;
      }
    }

    return false;
  }

  Hypergraph hg;
  Context context;
  std::unique_ptr<KMinusOneGainManager> manager;
};

class KMinusOneGainManagerNonParamTest : public KMinusOneGainManagerBaseTest, public Test {

};

class KMinusOneGainManagerParamTest : public KMinusOneGainManagerBaseTest, public TestWithParam<const char*> {
 protected:
  void SetUp() override {
    std::string graph_filename;
    PartitionID k;
    std::tie(graph_filename, k) = parseParams();
    loadGraph(graph_filename, k);
  }

  void ASSERT_THAT_ALL_GAINS_ARE_VALID() const {
    for (const HypernodeID& hn : hg.nodes()) {
      ASSERT_THAT_GAINS_ARE_VALID_FOR_HN(hn);
    }
  }

  void ASSERT_THAT_GAINS_ARE_VALID_FOR_HN(const HypernodeID hn) const {
    SCOPED_TRACE("bad gain values");
    std::vector<bool> checked_part(context.partition.k);
    checked_part[hg.partID(hn)] = true;

    if (!manager->isInternal(hn)) {
      for (const PartitionID& adjacent_part : manager->adjacentParts(hn)) {
        ASSERT_FALSE(checked_part[adjacent_part]);
        checked_part[adjacent_part] = true;
        ASSERT_TRUE(hypernodeIsConnectedToPart(hn, adjacent_part));
        ASSERT_THAT(manager->adjacentGain(hn, adjacent_part), Eq(gainInducedByHypergraph(hn, adjacent_part)))
                  << "bad gain adjacent gain value for " << V(hn) << ", " V(adjacent_part);
        ASSERT_THAT(manager->adjacentGain(hn, adjacent_part), Eq(manager->gain(hn, adjacent_part)));
      }
    }

    for (PartitionID nonadjacent_part = 0; nonadjacent_part < context.partition.k; ++nonadjacent_part) {
      if (checked_part[nonadjacent_part]) {
        continue;
      }

      ASSERT_FALSE(hypernodeIsConnectedToPart(hn, nonadjacent_part));
      ASSERT_THAT(manager->nonadjacentGain(hn), Eq(gainInducedByHypergraph(hn, nonadjacent_part)));
      ASSERT_THAT(manager->nonadjacentGain(hn), Eq(manager->gain(hn, nonadjacent_part)));
    }
  }

 private:
  std::pair<std::string, PartitionID> parseParams() const {
    const std::string param = GetParam();
    std::string graph_filename;
    PartitionID k;
    std::stringstream(param) >> graph_filename >> k;
    return {graph_filename, k};
  }
};

TEST_P(KMinusOneGainManagerParamTest, InitialGainsAreCorrect) {
  ASSERT_THAT_ALL_GAINS_ARE_VALID();
}

TEST_P(KMinusOneGainManagerParamTest, GainUpdateOfMovedHypernodeIsCorrect) {
  static constexpr std::size_t LIMIT = 128; // test at most this many hypernodes
  std::size_t num_hns_tested = 0;

  ASSERT_THAT_ALL_GAINS_ARE_VALID();

  for (const HypernodeID& hn : hg.nodes()) {
    const PartitionID from_part = hg.partID(hn);

    for (PartitionID to_part = 0; to_part < context.partition.k; ++to_part) {
      if (to_part == from_part) {
        continue;
      }

      hg.changeNodePart(hn, from_part, to_part);
      manager->updateAfterMovement(hn, from_part, to_part);
      ASSERT_THAT_GAINS_ARE_VALID_FOR_HN(hn);

      hg.changeNodePart(hn, to_part, from_part);
      manager->updateAfterMovement(hn, to_part, from_part);
      ASSERT_THAT_GAINS_ARE_VALID_FOR_HN(hn);
    }

    ++num_hns_tested;
    if (num_hns_tested == LIMIT) {
      break;
    }
  }
}

TEST_P(KMinusOneGainManagerParamTest, GainUpdateOfNeighborsOfMovedHypernodeIsCorrect) {
  static constexpr std::size_t LIMIT = 128; // move at most this many hypernodes
  std::size_t num_hns_tested = 0;
  ds::FastResetArray<bool> checked_neighbors(hg.initialNumNodes(), false);

  for (const HypernodeID& hn : hg.nodes()) {
    const PartitionID from_part = hg.partID(hn);

    for (PartitionID to_part = 0; to_part < context.partition.k; ++to_part) {
      if (to_part == from_part) {
        continue;
      }

      hg.changeNodePart(hn, from_part, to_part);
      manager->updateAfterMovement(hn, from_part, to_part);

      for (const HyperedgeID& he : hg.incidentEdges(hn)) {
        for (const HypernodeID& pin : hg.pins(he)) {
          if (checked_neighbors.get(pin)) {
            continue;
          }

          ASSERT_THAT_GAINS_ARE_VALID_FOR_HN(pin);
          checked_neighbors.set(pin, true);
        }
      }
      checked_neighbors.resetUsedEntries();

      hg.changeNodePart(hn, to_part, from_part);
      manager->updateAfterMovement(hn, to_part, from_part);

      for (const HyperedgeID& he : hg.incidentEdges(hn)) {
        for (const HypernodeID& pin : hg.pins(he)) {
          if (checked_neighbors.get(pin)) {
            continue;
          }

          ASSERT_THAT_GAINS_ARE_VALID_FOR_HN(pin);
          checked_neighbors.set(pin, true);
        }
      }
      checked_neighbors.resetUsedEntries();
    }

    ++num_hns_tested;
    if (num_hns_tested == LIMIT) {
      break;
    }
  }
}

TEST_P(KMinusOneGainManagerParamTest, GainUpdateOfNeighborsOfMovedHypernodeIsCorrectAfterMultipleMoves) {
  static constexpr std::size_t LIMIT = 128; // move at most this many hypernodes
  std::size_t num_hns_tested = 0;
  ds::FastResetArray<bool> checked_neighbors(hg.initialNumNodes(), false);

  for (const HypernodeID& hn : hg.nodes()) {
    PartitionID from_part = hg.partID(hn);

    for (PartitionID to_part = 0; to_part < context.partition.k; ++to_part) {
      if (to_part == from_part) {
        continue;
      }

      hg.changeNodePart(hn, from_part, to_part);
      manager->updateAfterMovement(hn, from_part, to_part);

      for (const HyperedgeID& he : hg.incidentEdges(hn)) {
        for (const HypernodeID& pin : hg.pins(he)) {
          if (checked_neighbors.get(pin)) {
            continue;
          }
          ASSERT_THAT_GAINS_ARE_VALID_FOR_HN(pin);
          checked_neighbors.set(pin, true);
        }
      }
      checked_neighbors.resetUsedEntries();

      from_part = to_part;
    }

    ++num_hns_tested;
    if (num_hns_tested == LIMIT) {
      break;
    }

    ASSERT_THAT_ALL_GAINS_ARE_VALID();
  }
}

TEST_P(KMinusOneGainManagerParamTest, GainUpdatesAreCorrectAfterRandomMovements) {
  static constexpr std::size_t LIMIT = 10'000; // number of moves
  auto& rand = Randomize::instance();

  for (std::size_t i = 0; i < LIMIT; ++i) {
    const HypernodeID hn = rand.getRandomInt(0, hg.initialNumNodes() - 1);
    const PartitionID from_part = hg.partID(hn);
    PartitionID to_part;
    do {
      to_part = rand.getRandomInt(0, context.partition.k - 1);
    } while (from_part == to_part);

    hg.changeNodePart(hn, from_part, to_part);
    manager->updateAfterMovement(hn, from_part, to_part);
  }

  ASSERT_THAT_ALL_GAINS_ARE_VALID();
}

class PQCallback {
 public:
  virtual void onConnectivityIncrease(HypernodeID hn, PartitionID part, Gain gain) = 0;
  virtual void onUpdate(HypernodeID hn, PartitionID part, Gain delta) = 0;
  virtual void onConnectivityDecrease(HypernodeID hn, PartitionID part, Gain gain) = 0;
};

class PQCallbackMock : public PQCallback {
 public:
  MOCK_METHOD3(onConnectivityIncrease, void(HypernodeID hn, PartitionID part, Gain gain));
  MOCK_METHOD3(onUpdate, void(HypernodeID hn, PartitionID part, Gain delta));
  MOCK_METHOD3(onConnectivityDecrease, void(HypernodeID hn, PartitionID part, Gain gain));
};

TEST_P(KMinusOneGainManagerParamTest, PQIsConsistentThroughCallbacks) {
  std::vector<std::vector<Gain>> pq(hg.initialNumNodes(), std::vector<Gain>(context.partition.k));
  constexpr Gain invalid_gain = std::numeric_limits<Gain>::max();
  constexpr std::size_t LIMIT = 10'000;

  for (const HypernodeID& hn : hg.nodes()) {
    for (PartitionID part = 0; part < context.partition.k; ++part) {
      if (hg.partID(hn) == part) {
        pq[hn][part] = invalid_gain;
      } else {
        pq[hn][part] = manager->gain(hn, part);
      }
    }
  }

  auto connectivity_inc_cb = [&pq](const auto& hn, const auto& part, const auto& gain) {
    pq[hn][part] = gain;
  };
  auto update_cb = [&pq, this](const auto& hn, const auto& part, const auto& delta) {
    if (part == Hypergraph::kInvalidPartition) { // nonadjacent gain changed
      for (PartitionID nonadjacent_part = 0; nonadjacent_part < context.partition.k; ++nonadjacent_part) {
        if (nonadjacent_part == hg.partID(hn)) {
          continue;
        }
        if (!hypernodeIsConnectedToPart(hn, nonadjacent_part) &&
            pq.at(hn).at(nonadjacent_part) != invalid_gain) { // this might be invalid_gain when updating moved_hn
          pq[hn][nonadjacent_part] += delta;
        }
      }
    } else { // gain for move to part changed
      ASSERT_THAT(pq.at(hn).at(part), Ne(invalid_gain));
      pq[hn][part] += delta;
    }
  };
  auto connectivity_dec_cb = [&pq](const auto& hn, const auto& part, const auto& gain) {
    ASSERT_THAT(pq.at(hn).at(part), Ne(invalid_gain));
    pq[hn][part] = gain;
  };
  auto move_cb = [&pq, this](const auto& hn, const auto& from_part, const auto& to_part) {
    pq[hn][from_part] = manager->gain(hn, from_part);
    pq[hn][to_part] = invalid_gain;
  };

  performRandomMoves(LIMIT, connectivity_inc_cb, update_cb, connectivity_dec_cb, move_cb);

  // validate
  for (const HypernodeID& hn : hg.nodes()) {
    for (PartitionID part = 0; part < context.partition.k; ++part) {
      if (hg.partID(hn) == part) {
        ASSERT_THAT(pq.at(hn).at(part), Eq(invalid_gain));
      } else {
        ASSERT_THAT(pq.at(hn).at(part), Eq(manager->gain(hn, part)));
      }
    }
  }
}

TEST_P(KMinusOneGainManagerParamTest, RandomCoarseningKeepsGainsCorrect) {
  constexpr std::size_t CONTRACTION_ITERATIONS = 3;
  constexpr std::size_t CONTRACTION_LIMIT = 1'000;
  std::vector<Hypergraph::Memento> contractions;
  auto& rand = Randomize::instance();
  std::vector<bool> matched(hg.initialNumNodes());
  std::vector<HypernodeID> current_hns;
  for (const HypernodeID& hn : hg.nodes()) {
    current_hns.push_back(hn);
  }

  std::size_t num_contraction = 0;
  for (std::size_t it = 0; it < CONTRACTION_ITERATIONS; ++it) {
    for (const HypernodeID& hn : current_hns) {
      if (matched[hn]) {
        continue;
      }

      for (const HyperedgeID& he : hg.incidentEdges(hn)) {
        for (const HypernodeID& pin : hg.pins(he)) {
          if (!matched[pin] && pin != hn && hg.partID(pin) == hg.partID(hn)) {
            matched[hn] = true;
            matched[pin] = true;
            contractions.push_back(hg.contract(hn, pin));
            break;
          }
        }
        if (matched[hn]) {
          break;
        }
      }
    }

    ++num_contraction;
    if (num_contraction == CONTRACTION_LIMIT) {
      break;
    }

    current_hns.clear();
    for (const HypernodeID& hn : hg.nodes()) {
      current_hns.push_back(hn);
    }
    matched.clear();
    matched.resize(hg.initialNumNodes(), false);
  }

  hg.initializeNumCutHyperedges();

  // init a new manager for the coarsest graph
  manager = std::make_unique<KMinusOneGainManager>(hg, context);
  manager->initialize();

  // finally, uncontract and make sure that everything stays valid
  ASSERT_THAT_ALL_GAINS_ARE_VALID();
  while (!contractions.empty()) {
    const auto& memento = contractions.back();

    manager->preUncontraction(memento.u);
    hg.uncontract(memento);
    manager->postUncontraction(memento.u, {memento.v});
    contractions.pop_back();
    ASSERT_THAT_ALL_GAINS_ARE_VALID();
  }
}

TEST_P(KMinusOneGainManagerParamTest, RollbackWorksWithRandomMovements) {
  constexpr std::size_t NUM_MOVES = 1000;

  std::vector<PartitionID> parts(hg.initialNumNodes());
  std::vector<std::vector<bool>> adjacent_to(hg.initialNumNodes(), std::vector<bool>(context.partition.k));
  std::vector<std::vector<Gain>> gain_to(hg.initialNumNodes(), std::vector<Gain>(context.partition.k));
  for (const HypernodeID& hn : hg.nodes()) {
    parts[hn] = hg.partID(hn);
    for (PartitionID part = 0; part < context.partition.k; ++part) {
      if (part == hg.partID(hn)) {
        continue;
      }

      adjacent_to[hn][part] = manager->isAdjacentTo(hn, part);
      gain_to[hn][part] = manager->gain(hn, part);
    }
  }

  manager->resetDelta();
  const auto& moves = performRandomMoves(NUM_MOVES);
  manager->rollbackDelta();

  // revert movements
  for (auto rit = moves.crbegin(); rit != moves.crend(); ++rit) {
    hg.changeNodePart(rit->hn, rit->to, rit->from);
  }

  ASSERT_THAT_ALL_GAINS_ARE_VALID();
  for (const HypernodeID& hn : hg.nodes()) {
    ASSERT_THAT(parts[hn], Eq(hg.partID(hn)));
    for (PartitionID part = 0; part < context.partition.k; ++part) {
      if (part == hg.partID(hn)) {
        continue;
      }

      ASSERT_THAT(adjacent_to[hn][part], Eq(manager->isAdjacentTo(hn, part)));
      ASSERT_THAT(gain_to[hn][part], Eq(manager->gain(hn, part)));
    }
  }
}

TEST_F(KMinusOneGainManagerNonParamTest, C17_CallbacksAreCorrectInArtificialTest) {
  loadGraph("test_instances/c17.hgr", 4);

  PQCallbackMock callback;
  auto connectivity_inc_cb = [&callback](const auto& hn, const auto& part, const auto& gain) {
    callback.onConnectivityIncrease(hn, part, gain);
  };
  auto update_cb = [&callback](const auto& hn, const auto& part, const auto& delta) {
    callback.onUpdate(hn, part, delta);
  };
  auto connectivity_dec_cb = [&callback](const auto& hn, const auto& part, const auto& gain) {
    callback.onConnectivityDecrease(hn, part, gain);
  };

  {
    EXPECT_CALL(callback, onConnectivityIncrease(0, 0, _));
    EXPECT_CALL(callback, onConnectivityIncrease(7, 0, _));
    EXPECT_CALL(callback, onConnectivityDecrease(1, 1, _));
    EXPECT_CALL(callback, onConnectivityDecrease(8, 1, _));

    constexpr HypernodeID hn = 2;
    constexpr PartitionID from = 1;
    constexpr PartitionID to = 0;
    hg.changeNodePart(hn, from, to);
    manager->updateAfterMovement(hn, from, to, connectivity_inc_cb, update_cb, connectivity_dec_cb);
  }

  // TODO extend this a bit if there are any errors with this
}

INSTANTIATE_TEST_CASE_P(GRAPH_C17_K_4, KMinusOneGainManagerParamTest, Values("test_instances/c17.hgr 4"));

INSTANTIATE_TEST_CASE_P(GRAPH_C17_K_8, KMinusOneGainManagerParamTest, Values("test_instances/c17.hgr 8"));

INSTANTIATE_TEST_CASE_P(GRAPH_C17_K_16, KMinusOneGainManagerParamTest, Values("test_instances/c17.hgr 16"));

INSTANTIATE_TEST_CASE_P(GRAPH_C3540_K_4, KMinusOneGainManagerParamTest, Values("test_instances/c3540.hgr 4"));

INSTANTIATE_TEST_CASE_P(GRAPH_C3540_K_8, KMinusOneGainManagerParamTest, Values("test_instances/c3540.hgr 8"));

INSTANTIATE_TEST_CASE_P(GRAPH_C3540_K_16, KMinusOneGainManagerParamTest, Values("test_instances/c3540.hgr 16"));

INSTANTIATE_TEST_CASE_P(GRAPH_C7552_K_16, KMinusOneGainManagerParamTest, Values("test_instances/c7552.hgr 16"));

INSTANTIATE_TEST_CASE_P(GRAPH_C7552_K_32, KMinusOneGainManagerParamTest, Values("test_instances/c7552.hgr 32"));
} // namespace dag
} // namespace kahypar