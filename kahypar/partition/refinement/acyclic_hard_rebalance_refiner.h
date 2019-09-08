#pragma once

#include <array>

#include "kahypar/datastructure/fast_reset_array.h"
#include "kahypar/datastructure/fast_reset_flag_array.h"
#include "kahypar/datastructure/sparse_map.h"
#include "kahypar/dag/quotient_graph.h"
#include "kahypar/definitions.h"
#include "kahypar/meta/template_parameter_to_string.h"
#include "kahypar/partition/context.h"
#include "kahypar/partition/metrics.h"
#include "kahypar/partition/refinement/i_refiner.h"
#include "kahypar/partition/refinement/km1_gain_manager.h"
#include "kahypar/datastructure/kway_priority_queue.h"
#include "kahypar/utils/randomize.h"
#include "kahypar/utils/timer.h"

namespace kahypar {
class AcyclicHardRebalanceRefiner final : public IRefiner {
 private:
  using KWayRefinementPQ = ds::KWayPriorityQueue<HypernodeID, Gain, std::numeric_limits<Gain>>;

  static constexpr bool debug = false;
  static constexpr HypernodeID hn_to_debug = 0;

  static constexpr std::size_t PREV = 0;
  static constexpr std::size_t NEXT = 1;
  static constexpr std::size_t TARGETS = 2;

  static constexpr std::size_t ACTIVATE = 1;
  static constexpr std::size_t DEACTIVATE = 2;

  static constexpr std::array<std::size_t, 2> kDirections{PREV, NEXT};

 public:
  AcyclicHardRebalanceRefiner(Hypergraph& hypergraph, const Context& context,
                              AdjacencyMatrixQuotientGraph<DFSCycleDetector>& qg,
                              KMinusOneGainManager& gain_manager) :
    _hg(hypergraph),
    _context(context),
    _qg(qg),
    _pq{KWayRefinementPQ(context.partition.k), KWayRefinementPQ(context.partition.k)},
    _fixtures{std::vector<HypernodeID>(), std::vector<HypernodeID>()},
    _state_changes{ds::FastResetArray<std::size_t>(hypergraph.initialNumNodes(), 0),
                   ds::FastResetArray<std::size_t>(hypergraph.initialNumNodes(), 0)},
    _updated_neighbors(_hg.initialNumNodes(), false),
    _gain_manager(gain_manager),
    _moved_from(_hg.initialNumNodes()) {}

  ~AcyclicHardRebalanceRefiner() override = default;

  AcyclicHardRebalanceRefiner(const AcyclicHardRebalanceRefiner&) = delete;
  AcyclicHardRebalanceRefiner& operator=(const AcyclicHardRebalanceRefiner&) = delete;

  AcyclicHardRebalanceRefiner(AcyclicHardRebalanceRefiner&&) = delete;
  AcyclicHardRebalanceRefiner& operator=(AcyclicHardRebalanceRefiner&&) = delete;

  const std::vector<Move>& moves() {
    return _moves;
  }

  double improvedImbalance() const {
    return _improved_imbalance;
  }

  void preUncontraction(const HypernodeID representant) override {
    if (_qg_changed) { // will refresh fixtures anyways
      return;
    }

    removeFixturesForHypernode(representant);
  }

  void postUncontraction(const HypernodeID representant, const std::vector<HypernodeID>&& partners) override {
    if (_qg_changed) { // will refresh fixtures anyways
      return;
    }

    ASSERT(partners.size() == 1, "Currently only supports pair contractions.");
    const HypernodeID partner = partners.front();

    addFixturesForHypernode(representant, partner);
    addFixturesForHypernode(partner, representant);
    processFixtureStateChanges();

    ASSERT([&]() {
      ASSERT_THAT_FIXTURES_ARE_CORRECT();
      ASSERT_THAT_HYPERNODES_CONTAINED_IN_PQS_ARE_CORRECT();
      return true;
    }());
  }

  void printSummary() const override {
    LOG << "[HardRebalance] Iterations with refresh:" << _num_refreshes;
    LOG << "[HardRebalance] Iterations without refresh:" << _num_no_refreshes;
    LOG << "[HardRebalance] Number of moves:" << _num_moves;
    LOG << "[HardRebalance] Number of moves in the last iteration:" << _num_moves_in_last_iteration;
    LOG << "[HardRebalance] Positive gain moves:" << _num_positive_gain_moves;
    LOG << "[HardRebalance] Zero gain moves:" << _num_zero_gain_moves;
    LOG << "[HardRebalance] Negative gain moves:" << _num_negative_gain_moves;
    LOG << "[HardRebalance] Improved imbalance by:" << _improved_imbalance;
    LOG << "[HardRebalance] Improved KM1 by:" << _improved_km1;
    LOG << "[HardRebalance] Number of times no path was found:" << _num_no_path;
    LOG << "[HardRebalance] Number of times a path was found:" << _num_found_path;
    LOG << "[HardRebalance] Number of times rebalance was successful:" << _num_success;
    LOG << "[HardRebalance] Perform moves time:" << _perform_moves_time;
    LOG << "[HardRebalance] Refine time:" << _refine_time;
    LOG << "[HardRebalance] Gains time:" << _gain_time;
    LOG << "[HardRebalance] Path time:" << _path_time;
  }

 private:
  void initializeImpl(const HyperedgeWeight max_gain) final {
    if (!_is_initialized) {
      _pq[PREV].initialize(_hg.initialNumNodes());
      _pq[NEXT].initialize(_hg.initialNumNodes());
      _is_initialized = true;
    }
    _hg.resetHypernodeState();
    refreshTopologicalOrdering();

    _num_refreshes = 0;
    _num_no_refreshes = 0;
    _num_moves = 0;
    _num_moves_in_last_iteration = 0;
    _num_positive_gain_moves = 0;
    _num_negative_gain_moves = 0;
    _num_zero_gain_moves = 0;
    _improved_imbalance = 0.0;
    _improved_km1 = 0;
    _num_no_path = 0;
    _num_found_path = 0;
    _num_success = 0;
    _perform_moves_time = 0.0;
    _refine_time = 0.0;
    _gain_time = 0.0;
    _path_time = 0.0;
  }

  void performMovesAndUpdateCacheImpl(const std::vector<Move>& moves,
                                      std::vector<HypernodeID>& refinement_nodes,
                                      const UncontractionGainChanges& changes) final {
    const auto& ordering = _qg.topologicalOrdering();
    if (_qg.changed() || _qg_changed) { // cached topological ordering changed, must reset PQs anyways
      _qg_changed = true;
      return;
    }

    const HighResClockTimepoint start_perform_moves = std::chrono::high_resolution_clock::now();

    _hg.resetHypernodeState();

    for (std::size_t i = moves.size(); i > 0; --i) {
      const Move& move = moves[i - 1];
      _hg.changeNodePart(move.hn, move.to, move.from);
    }

    for (const Move& move : moves) {
      for (const std::size_t& direction : kDirections) {
        if (_pq[direction].contains(move.hn, move.from)) {
          _pq[direction].remove(move.hn, move.from);
        }
      }
      _hg.changeNodePart(move.hn, move.from, move.to);
      this->move(move.hn, move.from, move.to);
    }

    const HighResClockTimepoint end_perform_moves = std::chrono::high_resolution_clock::now();
    _perform_moves_time += std::chrono::duration<double>(end_perform_moves - start_perform_moves).count();;

    ASSERT([&]() {
      ASSERT_THAT_FIXTURES_ARE_CORRECT();
      ASSERT_THAT_HYPERNODES_CONTAINED_IN_PQS_ARE_CORRECT();
      return true;
    }());
  }

  bool refineImpl(std::vector<HypernodeID>& refinement_nodes,
                  const std::array<HypernodeWeight, 2>&,
                  const UncontractionGainChanges&,
                  Metrics& best_metrics) final {
    const HighResClockTimepoint start_refine = std::chrono::high_resolution_clock::now();

    ASSERT(best_metrics.imbalance == metrics::imbalance(_hg, _context));
    ASSERT(best_metrics.km1 == metrics::km1(_hg));
    _hg.resetHypernodeState();
    _moves.clear();

    const auto& ordering = _qg.topologicalOrdering();
    if (_qg_changed || _qg.changed()) {
      refreshTopologicalOrdering();
      ++_num_refreshes;
    } else {
      ++_num_no_refreshes;
    }

    const double initial_imbalance = best_metrics.imbalance;
    const HyperedgeWeight initial_km1 = best_metrics.km1;

    ASSERT([&]() {
      ASSERT_THAT_FIXTURES_ARE_CORRECT();
      ASSERT_THAT_HYPERNODES_CONTAINED_IN_PQS_ARE_CORRECT();
      return true;
    }());
    ASSERT(_qg.isAcyclic());
    ASSERT(AdjacencyMatrixQuotientGraph<DFSCycleDetector>(_hg, _context).isAcyclic());

    bool stop_after_iteration = false;

    while (!stop_after_iteration && best_metrics.imbalance > _context.partition.epsilon) {
      const auto& path = selectPath();
      if (path.first == Hypergraph::kInvalidPartition || path.second == Hypergraph::kInvalidPartition) {
        ++_num_no_path;
        break;
      } else {
        ++_num_found_path;
      }

      const PartitionID& start = path.first;
      const PartitionID& end = path.second;
      const std::size_t direction = ordering[start] > ordering[end] ? PREV : NEXT;
      HyperedgeWeight km1_change = 0;

      for (PartitionID part = start; part != end; part = adjacentPart(part, direction)) {
        const PartitionID to_part = part;
        const PartitionID from_part = adjacentPart(to_part, direction);
        ASSERT(isAdjacent(from_part, to_part), V(from_part) << V(to_part) << V(ordering));
        ASSERT(_pq[reverse(direction)].isEnabled(from_part), V(from_part) << V(reverse(direction)));
        ASSERT(!_pq[reverse(direction)].empty(from_part), V(from_part) << V(reverse(direction)));
        HypernodeID max_gain_hn = 0;
        Gain max_gain = 0;
        _pq[reverse(direction)].deleteMaxFromPartition(max_gain_hn, max_gain, from_part);
        if (_pq[direction].contains(max_gain_hn, from_part)) {
          _pq[direction].remove(max_gain_hn, from_part);
        }
        ASSERT(_hg.partID(max_gain_hn) == from_part, V(max_gain_hn) << V(_hg.partID(max_gain_hn)) << V(from_part));

        if (_moved_from[max_gain_hn] > _hg.initialNumNodes()) {
          //stop_after_iteration = true;
        }

        ++_num_moves;
        if (_hg.initialNumNodes() == _hg.currentNumNodes()) {
          ++_num_moves_in_last_iteration;
        }
        if (max_gain > 0) {
          ++_num_positive_gain_moves;
        } else if (max_gain < 0) {
          ++_num_negative_gain_moves;
        } else {
          ++_num_zero_gain_moves;
        }
        km1_change -= max_gain;

        DBG << "Move HN" << max_gain_hn << "from" << from_part << "to" << to_part;

        // update QuotientGraph
        HighResClockTimepoint start_tp = std::chrono::high_resolution_clock::now();
        const bool success = _qg.testAndUpdateBeforeMovement(max_gain_hn, from_part, to_part);
        HighResClockTimepoint end_tp = std::chrono::high_resolution_clock::now();
        Timer::instance().add(_context, Timepoint::cycle_detector, std::chrono::duration<double>(end_tp - start_tp).count());
        ASSERT(success, V(max_gain_hn) << V(from_part) << V(to_part));

        // update Hypergraph
#ifdef KAHYPAR_USE_ASSERTIONS
        const auto km1_before_move = metrics::km1(_hg);
        const auto cached_gain_before_move = _gain_manager.gain(max_gain_hn, to_part);
#endif // KAHYPAR_USE_ASSERTIONS
        _hg.changeNodePart(max_gain_hn, from_part, to_part);
        ASSERT(km1_before_move - metrics::km1(_hg) == max_gain, V(km1_before_move)
          << V(metrics::km1(_hg))
          << V(max_gain)
          << V(cached_gain_before_move));

        // update gains
        const HighResClockTimepoint start_update_gc = std::chrono::high_resolution_clock::now();
        _gain_manager.updateAfterMovement(max_gain_hn, from_part, to_part);
        const HighResClockTimepoint end_update_gc = std::chrono::high_resolution_clock::now();
        _gain_time += std::chrono::duration<double>(end_update_gc - start_update_gc).count();;

        move(max_gain_hn, from_part, to_part);
        _moves.emplace_back(max_gain_hn, from_part, to_part);
      }

      const double new_imbalance = metrics::imbalance(_hg, _context);
      best_metrics.imbalance = new_imbalance;
      best_metrics.km1 += km1_change;
      ASSERT(best_metrics.km1 == metrics::km1(_hg));
      _gain_manager.resetDelta();
    }

    ASSERT([&]() {
      ASSERT_THAT_FIXTURES_ARE_CORRECT();
      ASSERT_THAT_HYPERNODES_CONTAINED_IN_PQS_ARE_CORRECT();
      return true;
    }());
    ASSERT(_qg.isAcyclic());
    ASSERT(AdjacencyMatrixQuotientGraph<DFSCycleDetector>(_hg, _context).isAcyclic());
    // ASSERT(best_metrics.imbalance <= _context.partition.epsilon); -- might fail if the graph is cyclic

    DBG << "Imbalance after hard rebalance refiner:" << best_metrics.imbalance;
    ASSERT(best_metrics.km1 == metrics::km1(_hg));
    ASSERT(best_metrics.imbalance == metrics::imbalance(_hg, _context));
    _improved_imbalance += initial_imbalance - best_metrics.imbalance;
    _improved_km1 += initial_km1 - best_metrics.km1;

    if (best_metrics.imbalance <= _context.partition.epsilon) {
      ++_num_success;
    }

    for (const auto& hn : _moved_from_entries) {
      _moved_from[hn] = 0;
    }
    _moved_from_entries.clear();

    const HighResClockTimepoint end_refine = std::chrono::high_resolution_clock::now();
    _refine_time += std::chrono::duration<double>(end_refine - start_refine).count();;

    return false;
  }

  void logPQStats() const {
    DBG << "Topological ordering:" << _qg.topologicalOrdering();
    DBG << "Inverse topological ordering:" << _qg.inverseTopologicalOrdering();
    for (const auto& part : _qg.inverseTopologicalOrdering()) {
      DBG << "\tPartition" << part << "-->" << _pq[PREV].size(part) << "/" << _pq[NEXT].size(part) << "/"
          << _hg.partWeight(part) << "/" << _hg.partSize(part);
    }
  }

  /*********************************************************************************
   * SCHEDULING METHODS
   *********************************************************************************/
  void refreshTopologicalOrdering() {
    DBG << "Re-initialize PQs in AcyclicHardRebalanceRefiner";
    _fixtures[PREV].clear();
    _fixtures[NEXT].clear();
    _fixtures[PREV].resize(_hg.initialNumNodes());
    _fixtures[NEXT].resize(_hg.initialNumNodes());
    _pq[PREV].clear();
    _pq[NEXT].clear();

    initializeFixtures();
    ASSERT([&]() {
      ASSERT_THAT_FIXTURES_ARE_CORRECT();
      return true;
    }());

    processFixtureStateChanges();
    ASSERT([&]() {
      ASSERT_THAT_HYPERNODES_CONTAINED_IN_PQS_ARE_CORRECT();
      return true;
    }());

    _qg_changed = false;
  }

  std::size_t reverse(const std::size_t direction) const {
    return 1 - direction;
  }

  bool isOverloaded(const PartitionID part) const {
    return _hg.partWeight(part) > _context.partition.max_part_weights[part];
  }

  bool isUnderloaded(const PartitionID part) const {
    return _hg.partWeight(part) < (1.0 - _context.partition.epsilon) * _context.partition.perfect_balance_part_weights[part];
  }

  std::pair<PartitionID, PartitionID> selectPath() const {
    const HighResClockTimepoint start_path = std::chrono::high_resolution_clock::now();

    const auto& ordering = _qg.topologicalOrdering();
    const auto& inverse_ordering = _qg.inverseTopologicalOrdering();

    Gain current_gain = 0;
    Gain best_gain = std::numeric_limits<Gain>::min();
    PartitionID best_start = Hypergraph::kInvalidPartition;
    PartitionID best_end = Hypergraph::kInvalidPartition;
    std::size_t best_direction = 0;

    for (PartitionID from_id = 0; from_id < _context.partition.k; ++from_id) {
      const PartitionID from_part = inverse_ordering[from_id];
      if (!isOverloaded(from_part)) {
        continue;
      }

      current_gain = 0;

      for (PartitionID to_id_prime = from_id; to_id_prime > 0; --to_id_prime) {
        const PartitionID to_id = to_id_prime - 1;
        ASSERT(from_id != to_id);
        const PartitionID to_part = inverse_ordering[to_id];
        ASSERT(from_part != to_part);
        const PartitionID prev_part = inverse_ordering[to_id_prime];

        if (_pq[PREV].empty(prev_part)) { // path not feasible
          break;
        }
        const Gain gain = _pq[PREV].maxKey(prev_part);
        current_gain += gain;

        if (current_gain > best_gain && isUnderloaded(to_part)) {
          best_gain = current_gain;
          best_start = from_part;
          best_end = to_part;
          best_direction = PREV;
        }
      }

      current_gain = 0;

      for (PartitionID to_id_prime = from_id; to_id_prime + 1 < _context.partition.k; ++to_id_prime) {
        const PartitionID to_id = to_id_prime + 1;
        ASSERT(to_id != from_id);
        const PartitionID to_part = inverse_ordering[to_id];
        ASSERT(to_part != from_part);
        const PartitionID prev_part = inverse_ordering[to_id_prime];
        ASSERT(adjacentPart(prev_part, NEXT) == to_part);
        if (_pq[NEXT].empty(prev_part)) { // path not feasible
          break;
        }
        const Gain gain = _pq[NEXT].maxKey(prev_part);
        current_gain += gain;

        if (current_gain > best_gain && isUnderloaded(to_part)) {
          best_gain = current_gain;
          best_start = from_part;
          best_end = to_part;
          best_direction = NEXT;
        }
      }
    }

    DBG << "Path:" << best_start << "-->" << best_end << "/" << best_direction << "/" << best_gain;

    ASSERT([&]() {
      if (best_start != Hypergraph::kInvalidPartition || best_end != Hypergraph::kInvalidPartition) {
        ASSERT(best_start >= 0 && best_start < _context.partition.k);
        ASSERT(best_end >= 0 && best_end < _context.partition.k);
        ASSERT(best_start != best_end);
        if (best_direction == PREV) {
          ASSERT(ordering[best_start] > ordering[best_end]);
        } else {
          ASSERT(ordering[best_end] > ordering[best_start]);
        }

        const std::size_t end_to_start = reverse(best_direction);
        for (PartitionID part = best_end; part != best_start; part = adjacentPart(part, end_to_start)) {
          DBG << V(part) << V(best_end) << V(best_start) << V(best_direction) << V(adjacentPart(part, end_to_start)) << V(end_to_start);
          ASSERT(!_pq[best_direction].empty(adjacentPart(part, end_to_start)));
        }
      }
      return true;
    }());

    const HighResClockTimepoint end_path = std::chrono::high_resolution_clock::now();
    const_cast<AcyclicHardRebalanceRefiner*>(this)->_path_time += std::chrono::duration<double>(end_path - start_path).count();

    return {best_end, best_start};
  }

  /*********************************************************************************
   * GAIN UPDATE METHODS
   *********************************************************************************/
  void move(const HypernodeID moved_hn, const PartitionID from_part, const PartitionID to_part) {
    if (_moved_from[moved_hn] == 0) {
      _moved_from_entries.emplace_back(moved_hn);
    }
    _moved_from[moved_hn]++;

    updateFixtures(moved_hn, from_part, to_part);
    processFixtureStateChanges();

    const HighResClockTimepoint start_update_gains = std::chrono::high_resolution_clock::now();

    _updated_neighbors.resetUsedEntries();

    updateHypernodeGain(moved_hn);
    _updated_neighbors.set(moved_hn, true);

    for (const HyperedgeID& he : _hg.incidentEdges(moved_hn)) {
      for (const HypernodeID& pin : _hg.pins(he)) {
        if (!_updated_neighbors.get(pin)) {
          updateHypernodeGain(pin);
          _updated_neighbors.set(pin, true);
        }
      }
    }

    const HighResClockTimepoint end_update_gains = std::chrono::high_resolution_clock::now();
    _gain_time += std::chrono::duration<double>(end_update_gains - start_update_gains).count();;
  }

  void updateHypernodeGain(const HypernodeID hn) {
    DBGC(hn == hn_to_debug) << "updateHypernodeGain(" << hn << ")";

    for (const std::size_t direction : kDirections) {
      if (_fixtures[direction][hn] != 0) {
        continue;
      }

      const PartitionID from_part = _hg.partID(hn);
      const PartitionID to_part = adjacentPart(from_part, direction);
      if (to_part == Hypergraph::kInvalidPartition) {
        continue;
      }

      ASSERT(_pq[direction].contains(hn, from_part));
      DBGC(hn == hn_to_debug) << "Update gain from" << _pq[direction].key(hn, from_part)
                              << "to" << _gain_manager.gain(hn, to_part)
                              << "for HN" << hn << "and parts" << from_part << " / " << to_part;
      _pq[direction].updateKey(hn, from_part, _gain_manager.gain(hn, to_part));
    }
  }

  std::size_t directionToAdjacent(const PartitionID from, const PartitionID to) {
    return adjacentPart(from, PREV) == to ? PREV : NEXT;
  }

  void updateFixtures(const HypernodeID moved_hn, const PartitionID from_part, const PartitionID to_part) {
    DBGC(moved_hn == hn_to_debug) << "updateFixtures(" << V(moved_hn) << "," << V(from_part) << "," << V(to_part) << ")";
    _fixtures[PREV][moved_hn] = 0;
    _fixtures[NEXT][moved_hn] = 0;

    for (const HyperedgeID& he : _hg.incidentHeadEdges(moved_hn)) {
      for (const HypernodeID& tail : _hg.tails(he)) {
        if (from_part == _hg.partID(tail)) {
          ASSERT(_fixtures[NEXT][tail] > 0);
          --_fixtures[NEXT][tail];
          if (_fixtures[NEXT][tail] == 0) {
            _state_changes[NEXT].set(tail, ACTIVATE);
          }
        }
        if (to_part == _hg.partID(tail)) {
          ++_fixtures[NEXT][tail];
          ++_fixtures[PREV][moved_hn];
          if (_fixtures[NEXT][tail] == 1) {
            _state_changes[NEXT].set(tail, DEACTIVATE);
          }
        }
      }
    }
    for (const HyperedgeID& he : _hg.incidentTailEdges(moved_hn)) {
      for (const HypernodeID& head : _hg.heads(he)) {
        if (from_part == _hg.partID(head)) {
          ASSERT(_fixtures[PREV][head] > 0);
          --_fixtures[PREV][head];
          if (_fixtures[PREV][head] == 0) {
            _state_changes[PREV].set(head, ACTIVATE);
          }
        }
        if (to_part == _hg.partID(head)) {
          ++_fixtures[PREV][head];
          ++_fixtures[NEXT][moved_hn];
          if (_fixtures[PREV][head] == 1) {
            _state_changes[PREV].set(head, DEACTIVATE);
          }
        }
      }
    }

    if (_fixtures[PREV][moved_hn] == 0) {
      DBGC(moved_hn == hn_to_debug) << "_state_changes in PREV to ACTIVATE for" << V(moved_hn);
      _state_changes[PREV].set(moved_hn, ACTIVATE);
    }
    if (_fixtures[NEXT][moved_hn] == 0) {
      DBGC(moved_hn == hn_to_debug) << "_state_changes in NEXT to ACTIVATE for" << V(moved_hn);
      _state_changes[NEXT].set(moved_hn, ACTIVATE);
    }
  }

  bool isFirstPartition(const PartitionID part) const {
    return _qg.topologicalOrdering()[part] == 0;
  }

  bool isLastPartition(const PartitionID part) const {
    return _qg.topologicalOrdering()[part] == _context.partition.k - 1;
  }

  bool isAdjacent(const PartitionID a, const PartitionID b) const {
    const auto& ordering = _qg.topologicalOrdering();
    return std::abs(static_cast<int>(ordering[a]) - static_cast<int>(ordering[b])) == 1;
  }

  bool isValidDirection(const PartitionID part, std::size_t direction) const {
    return !((direction == NEXT && isLastPartition(part)) || (direction == PREV && isFirstPartition(part)));
  }

  void processFixtureStateChanges() {
    DBG << "processFixtureStateChanges(): " << V(_state_changes[PREV].size()) << V(_state_changes[NEXT].size());

    for (const std::size_t& direction : kDirections) {
      for (const auto& hn : _state_changes[direction].usedEntries()) {
        const PartitionID part = _hg.partID(hn);
        if (!isValidDirection(part, direction)) {
          DBG << "Invalid direction" << V(part) << V(direction);
          continue;
        }

        if (_state_changes[direction].get(hn) == ACTIVATE) {
          const PartitionID to_part = adjacentPart(part, direction);
          const Gain gain = _gain_manager.gain(hn, to_part);
          DBG << "Activate" << V(hn) << V(direction) << V(part) << V(to_part) << V(gain);
          if (_pq[direction].contains(hn, part)) {
            _pq[direction].updateKey(hn, part, gain);
          } else {
            _pq[direction].insert(hn, part, gain);
          }
          _pq[direction].enablePart(part);
        } else if (_state_changes[direction].get(hn) == DEACTIVATE) {
          DBG << "Deactivate HN" << V(hn) << V(direction) << V(part);
          if (_pq[direction].contains(hn, part)) {
            _pq[direction].remove(hn, part);
          }
          if (_pq[direction].empty(part)) {
            _pq[direction].disablePart(part);
          }
        }
      }
    }

    for (const auto& direction : kDirections) {
      _state_changes[direction].resetUsedEntries();
    }
  }

  void initializeFixtures() {
    for (const HypernodeID& hn : _hg.nodes()) {
      initializeFixturesFor(hn);
    }
  }

  void initializeFixturesFor(const HypernodeID hn) {
    const PartitionID part = _hg.partID(hn);
    for (const auto& direction : kDirections) {
      if (_fixtures[direction][hn] == 0 && _pq[direction].contains(hn, part)) {
        _state_changes[direction].set(hn, DEACTIVATE);
      }
      _fixtures[direction][hn] = 0;
    }

    for (const HyperedgeID& he : _hg.incidentHeadEdges(hn)) {
      for (const HypernodeID& tail : _hg.tails(he)) {
        if (_hg.partID(hn) == _hg.partID(tail)) {
          ++_fixtures[PREV][hn];
        }
      }
    }
    for (const HyperedgeID& he : _hg.incidentTailEdges(hn)) {
      for (const HypernodeID& head : _hg.heads(he)) {
        if (_hg.partID(hn) == _hg.partID(head)) {
          ++_fixtures[NEXT][hn];
        }
      }
    }

    DBGC(hn == hn_to_debug) << "Initial fixtures for" << V(hn) << V(_fixtures[PREV][hn]) << V(_fixtures[NEXT][hn]);

    if (_fixtures[PREV][hn] == 0) {
      DBGC(hn == hn_to_debug) << "_state_changes in PREV to ACTIVATE for" << V(hn);
      _state_changes[PREV].set(hn, ACTIVATE);
    }
    if (_fixtures[NEXT][hn] == 0) {
      DBGC(hn == hn_to_debug) << "_state_changes in NEXT to ACTIVATE for" << V(hn);
      _state_changes[NEXT].set(hn, ACTIVATE);
    }
  }

  void addFixturesForHypernode(const HypernodeID hn, const HypernodeID ignore) {
    for (const HyperedgeID& he : _hg.incidentTailEdges(hn)) {
      for (const HypernodeID& head : _hg.heads(he)) {
        if (_hg.partID(hn) == _hg.partID(head)) { // (hn --> head)
          if (head != ignore) {
            incrementFixtures(head, PREV);
          }
          incrementFixtures(hn, NEXT);
        }
      }
    }
    for (const HyperedgeID& he : _hg.incidentHeadEdges(hn)) {
      for (const HypernodeID& tail : _hg.tails(he)) {
        if (_hg.partID(hn) == _hg.partID(tail)) { // (tail -> hn)
          if (tail != ignore) {
            incrementFixtures(tail, NEXT);
          }
          incrementFixtures(hn, PREV);
        }
      }
    }

    if (_fixtures[PREV][hn] == 0) {
      DBGC(hn == hn_to_debug) << "State changed to ACTIVE for PREV /" << hn;
      _state_changes[PREV].set(hn, ACTIVATE);
    }
    if (_fixtures[NEXT][hn] == 0) {
      DBGC(hn == hn_to_debug) << "State changed to ACTIVE for NEXT /" << hn;
      _state_changes[NEXT].set(hn, ACTIVATE);
    }
  }

  void removeFixturesForHypernode(const HypernodeID hn) {
    for (const HyperedgeID& he : _hg.incidentTailEdges(hn)) {
      for (const HypernodeID& head : _hg.heads(he)) {
        if (_hg.partID(hn) == _hg.partID(head)) { // (hn --> head)
          decrementFixtures(head, PREV);
          decrementFixtures(hn, NEXT);
        }
      }
    }
    for (const HyperedgeID& he : _hg.incidentHeadEdges(hn)) {
      for (const HypernodeID& tail : _hg.tails(he)) {
        if (_hg.partID(hn) == _hg.partID(tail)) { // (tail -> hn)
          decrementFixtures(tail, NEXT);
          decrementFixtures(hn, PREV);
        }
      }
    }

    ASSERT(_fixtures[PREV][hn] == 0);
    ASSERT(_fixtures[NEXT][hn] == 0);
    _state_changes[PREV].set(hn, DEACTIVATE);
    _state_changes[NEXT].set(hn, DEACTIVATE);
  }

  void incrementFixtures(const HypernodeID hn, const std::size_t direction) {
    if (_fixtures[direction][hn] == 0) {
      _state_changes[direction].set(hn, DEACTIVATE);
    }
    ++_fixtures[direction][hn];
  }

  void decrementFixtures(const HypernodeID hn, const std::size_t direction, const std::size_t by = 1) {
    if (_fixtures[direction][hn] == 1) {
      _state_changes[direction].set(hn, ACTIVATE);
    }
    ASSERT(_fixtures[direction][hn] >= by);
    _fixtures[direction][hn] -= by;
  }

  PartitionID adjacentPart(const PartitionID part, const std::size_t direction) const {
    const auto& ordering = _qg.topologicalOrdering();
    const auto& inverse_ordering = _qg.inverseTopologicalOrdering();

    switch (direction) {
      case PREV:
        return ordering[part] == 0
               ? Hypergraph::kInvalidPartition
               : inverse_ordering[ordering[part] - 1];

      case NEXT:
        return ordering[part] == _context.partition.k - 1
               ? Hypergraph::kInvalidPartition
               : inverse_ordering[ordering[part] + 1];

      default:
        ASSERT(false, "Invalid direction:" << direction);
        return 0;
    }
  }

  bool hypernodeIsConnectedToPart(const HypernodeID pin, const PartitionID part) const {
    for (const HyperedgeID& he : _hg.incidentEdges(pin)) {
      if (_hg.pinCountInPart(he, part) > 0) {
        return true;
      }
    }
    return false;
  }

#ifdef KAHYPAR_USE_ASSERTIONS
  void ASSERT_THAT_FIXTURES_ARE_CORRECT() const {
    std::array<std::vector<HypernodeID>, 2> fixtures{std::vector<HypernodeID>(_hg.initialNumNodes()),
                                                     std::vector<HypernodeID>(_hg.initialNumNodes())};

    for (const HypernodeID& hn : _hg.nodes()) {
      for (const HyperedgeID& he : _hg.incidentHeadEdges(hn)) {
        for (const HypernodeID& tail : _hg.tails(he)) {
          if (_hg.partID(hn) == _hg.partID(tail)) {
            ++fixtures[PREV][hn];
          }
        }
      }
      for (const HyperedgeID& he : _hg.incidentTailEdges(hn)) {
        for (const HypernodeID& head : _hg.heads(he)) {
          if (_hg.partID(hn) == _hg.partID(head)) {
            ++fixtures[NEXT][hn];
          }
        }
      }
    }

    for (const HypernodeID& hn : _hg.nodes()) {
      for (const std::size_t& direction : kDirections) {
        ASSERT (fixtures[direction][hn] == _fixtures[direction][hn]);
      }
    }
  }

  void ASSERT_THAT_HYPERNODES_CONTAINED_IN_PQS_ARE_CORRECT() const {
    for (const HypernodeID& hn : _hg.nodes()) {
      for (const std::size_t& direction : kDirections) {
        if (_fixtures[direction][hn] == 0
            && (direction != PREV || !isFirstPartition(_hg.partID(hn)))
            && (direction != NEXT || !isLastPartition(_hg.partID(hn)))) {
          ASSERT(_pq[direction].contains(hn));
          ASSERT(_pq[direction].contains(hn, _hg.partID(hn)));

          const PartitionID from_part = _hg.partID(hn);
          const PartitionID to_part = adjacentPart(from_part, direction);
          ASSERT(to_part != Hypergraph::kInvalidPartition);
          ASSERT(_pq[direction].key(hn, from_part) == _gain_manager.gain(hn, to_part),
            "In PQ:" << _pq[direction].key(hn, _hg.partID(hn))
            << "Cached:" << _gain_manager.gain(hn, to_part)
            << "Real:" << gainInducedByHypergraph(hn, to_part)
            << "For:" << V(hn) << V(from_part) << V(to_part));
        } else {
          ASSERT(!_pq[direction].contains(hn));
          ASSERT(!_pq[direction].contains(hn, _hg.partID(hn)));
        }
      }
    }
  }

  Gain gainInducedByHypergraph(const HypernodeID hn, const PartitionID target_part) const {
    Gain gain = 0;
    for (const HyperedgeID& he : _hg.incidentEdges(hn)) {
      ASSERT(_hg.edgeSize(he) > 1, V(he));
      gain += gainInducedByHyperedge(hn, he, target_part);
    }
    return gain;
  }

  Gain gainInducedByHyperedge(const HypernodeID hn, const HyperedgeID he, const PartitionID target_part) const {
    const HypernodeID pins_in_source_part = _hg.pinCountInPart(he, _hg.partID(hn));
    const HypernodeID pins_in_target_part = _hg.pinCountInPart(he, target_part);
    const HyperedgeWeight he_weight = _hg.edgeWeight(he);
    Gain gain = pins_in_source_part == 1 ? he_weight : 0;
    gain -= pins_in_target_part == 0 ? he_weight : 0;
    return gain;
  }
#endif

  Hypergraph& _hg;
  const Context& _context;
  AdjacencyMatrixQuotientGraph<DFSCycleDetector>& _qg;
  std::array<KWayRefinementPQ, TARGETS> _pq;
  std::array<std::vector<HypernodeID>, TARGETS> _fixtures;
  std::array<ds::FastResetArray<std::size_t>, 2> _state_changes;
  std::vector<Move> _moves;
  std::vector<HypernodeID> _moved_hns;
  ds::FastResetArray<bool> _updated_neighbors;
  KMinusOneGainManager& _gain_manager;
  std::vector<std::size_t> _moved_from;
  std::vector<HypernodeID> _moved_from_entries{};

  bool _qg_changed{false};

  std::size_t _num_refreshes{0};
  std::size_t _num_no_refreshes{0};
  std::size_t _num_moves{0};
  std::size_t _num_positive_gain_moves{0};
  std::size_t _num_negative_gain_moves{0};
  std::size_t _num_zero_gain_moves{0};
  std::size_t _num_found_path{0};
  std::size_t _num_no_path{0};
  std::size_t _num_success{0};
  std::size_t _num_moves_in_last_iteration{0};
  HyperedgeWeight _improved_km1{0};
  double _improved_imbalance{0.0};
  double _perform_moves_time{0.0};
  double _refine_time{0.0};
  double _gain_time{0.0};
  double _path_time{0.0};
};

constexpr std::array<std::size_t, 2> AcyclicHardRebalanceRefiner::kDirections;
} // namespace kahypar