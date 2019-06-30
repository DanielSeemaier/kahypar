#pragma once

#include <array>

#include "kahypar/dag/quotient_graph.h"

namespace kahypar {
class AcyclicHardRebalanceRefiner final : public IRefiner {
 private:
  using GainCache = KwayGainCache<Gain>;
  using KWayRefinementPQ = ds::KWayPriorityQueue<HypernodeID, Gain, std::numeric_limits<Gain>>;

  static constexpr bool debug = false;
  static constexpr HypernodeID hn_to_debug = 1261;

  static constexpr std::size_t PREV = 0;
  static constexpr std::size_t NEXT = 1;
  static constexpr std::size_t TARGETS = 2;

  static constexpr std::size_t ACTIVATE = 1;
  static constexpr std::size_t DEACTIVATE = 2;

  static constexpr std::array<std::size_t, 2> kDirections{PREV, NEXT};

  struct PinState {
    bool one_pin_in_from_part_before;
    bool one_pin_in_to_part_after;
    bool two_pins_in_from_part_before;
    bool two_pins_in_to_part_after;

    PinState(const bool one_in_from_before, const bool one_in_to_after,
             const bool two_in_from_before, const bool two_in_to_after) :
      one_pin_in_from_part_before(one_in_from_before),
      one_pin_in_to_part_after(one_in_to_after),
      two_pins_in_from_part_before(two_in_from_before),
      two_pins_in_to_part_after(two_in_to_after) {}
  };

 public:
  AcyclicHardRebalanceRefiner(Hypergraph& hypergraph, const Context& context) :
    _hg(hypergraph),
    _context(context),
    _pq{KWayRefinementPQ(context.partition.k), KWayRefinementPQ(context.partition.k)},
    _fixtures{std::vector<HypernodeID>(), std::vector<HypernodeID>()},
    _state_changes{ds::FastResetArray<std::size_t>(hypergraph.initialNumNodes(), 0),
                   ds::FastResetArray<std::size_t>(hypergraph.initialNumNodes(), 0)},
    _gain_cache(_hg.initialNumNodes(), _context.partition.k),
    _tmp_gains(_context.partition.k, 0),
    _nonadjacent_gain_cache(_hg.initialNumNodes()),
    _new_adjacent_part(_hg.initialNumNodes(), Hypergraph::kInvalidPartition) {}

  ~AcyclicHardRebalanceRefiner() override = default;

  AcyclicHardRebalanceRefiner(const AcyclicHardRebalanceRefiner&) = delete;
  AcyclicHardRebalanceRefiner& operator=(const AcyclicHardRebalanceRefiner&) = delete;

  AcyclicHardRebalanceRefiner(AcyclicSoftRebalanceRefiner&&) = delete;
  AcyclicHardRebalanceRefiner& operator=(AcyclicHardRebalanceRefiner&&) = delete;

  void changeEpsilon(const double epsilon) {
    _context.partition.epsilon = epsilon;
    _context.setupPartWeights(_hg.totalWeight());
  }

  const std::vector<Move>& moves() {
    return _moves;
  }

  void enableRefinementNodes(std::vector<HypernodeID>& refinement_nodes) {
    LOG << "enableRefinementNodes(" << refinement_nodes << ")";
    _hg.resetHypernodeState();
    for (const HypernodeID& hn : refinement_nodes) {
      DBGC(hn == hn_to_debug) << "Reset" << V(hn);
      _gain_cache.clear(hn);
      initializeGainCacheFor(hn);

      for (const HyperedgeID& he : _hg.incidentEdges(hn)) {
        for (const HypernodeID& pin : _hg.pins(he)) {
          initializeFixturesFor(pin);
        }
      }
    }

    processFixtureStateChanges();

    ASSERT([&]() {
      ASSERT_THAT_FIXTURES_ARE_CONSISTENT();
      return true;
    }());
  }

 private:
  void initializeImpl(const HyperedgeWeight max_gain) final {
    LOG << "initializeImpl()";
    if (!_is_initialized) {
      _pq[PREV].initialize(_hg.initialNumNodes());
      _pq[NEXT].initialize(_hg.initialNumNodes());
      _is_initialized = true;
    }
    _hg.resetHypernodeState();
    _gain_cache.clear();
    initializeGainCache();
    refreshTopologicalOrdering();
  }

  void performMovesAndUpdateCacheImpl(const std::vector<Move>& moves,
                                      std::vector<HypernodeID>& refinement_nodes,
                                      const UncontractionGainChanges& changes) final {
    LOG << "performMovesAndUpdateCacheImpl()";
    _hg.resetHypernodeState();

    for (std::size_t i = moves.size(); i > 0; --i) {
      const Move& move = moves[i - 1];
      _hg.changeNodePart(move.hn, move.to, move.from);
    }

    for (const Move& move : moves) {
      DBG << "move(" << V(move.hn) << V(move.from) << V(move.to) << ")";
      if (!_hg.active(move.hn)) {
        _hg.activate(move.hn);
      }
      for (const std::size_t& direction : kDirections) {
        if (_pq[direction].contains(move.hn, move.from)) {
          _pq[direction].remove(move.hn, move.from);
        }
      }
      _hg.changeNodePart(move.hn, move.from, move.to);
      this->move(move.hn, move.from, move.to);
    }

    ASSERT([&]() {
      ASSERT_THAT_FIXTURES_ARE_CONSISTENT();
      return true;
    }());
  }

  bool refineImpl(std::vector<HypernodeID>& refinement_nodes,
                  const std::array<HypernodeWeight, 2>&,
                  const UncontractionGainChanges&,
                  Metrics& best_metrics) final {
    LOG << "refineImpl()";
    LOG << "refineImpl(" << refinement_nodes << "): imbalance" << best_metrics.imbalance;
    ASSERT(best_metrics.imbalance == metrics::imbalance(_hg, _context));
    _moves.clear();
    _moved_hns.clear();
    _hg.resetHypernodeState();

    if (!orderingStillTopological()) {
      refreshTopologicalOrdering();
    } else {
      activateMovableHNs();
    }

    while (best_metrics.imbalance > _context.partition.epsilon) {
      DBG << "while(" << V(best_metrics.imbalance) << V(_context.partition.epsilon) << ")";
      ASSERT([&]() {
        ASSERT_THAT_FIXTURES_ARE_CONSISTENT();
        ASSERT_THAT_GAIN_CACHE_IS_VALID();
        return true;
      }());
      logPQStats();
      const auto& path = selectPath();
      const PartitionID& start = path.first;
      const PartitionID& end = path.second;
      const std::size_t direction = _ordering[start] > _ordering[end] ? PREV : NEXT;

      for (PartitionID part = start; part != end; part = adjacentPart(part, direction)) {
        const PartitionID to_part = part;
        const PartitionID from_part = adjacentPart(to_part, direction);
        ASSERT(isAdjacent(from_part, to_part), V(from_part) << V(to_part) << V(_ordering) << V(_inverse_ordering));
        ASSERT(_pq[reverse(direction)].isEnabled(from_part), V(from_part) << V(reverse(direction)));
        ASSERT(!_pq[reverse(direction)].empty(from_part), V(from_part) << V(reverse(direction)));
        HypernodeID max_gain_hn = 0;
        Gain max_gain = 0;
        _pq[reverse(direction)].deleteMaxFromPartition(max_gain_hn, max_gain, from_part);
        if (_pq[direction].contains(max_gain_hn, from_part)) {
          _pq[direction].remove(max_gain_hn, from_part);
        }
        ASSERT(_hg.active(max_gain_hn), V(max_gain_hn) << "not active");
        ASSERT(_hg.partID(max_gain_hn) == from_part, V(max_gain_hn) << V(_hg.partID(max_gain_hn)) << V(from_part));

        DBG << "refine.while(move:" << V(max_gain_hn) << V(max_gain) << ")";
        ASSERT([&]() {
          ASSERT_THAT_FIXTURES_ARE_CONSISTENT();
          return true;
        }());

        _hg.changeNodePart(max_gain_hn, from_part, to_part);
        move(max_gain_hn, from_part, to_part);
        _moves.emplace_back(max_gain_hn, from_part, to_part);
        _moved_hns.push_back(max_gain_hn);
      }

      const double new_imbalance = metrics::imbalance(_hg, _context);
      best_metrics.imbalance = new_imbalance;
    }

    ASSERT([&]() {
      ASSERT_THAT_FIXTURES_ARE_CONSISTENT();
      ASSERT_THAT_GAIN_CACHE_IS_VALID();
      return true;
    }());
    ASSERT(best_metrics.imbalance <= _context.partition.epsilon);

    DBG << "Imbalance after refineImpl():" << best_metrics.imbalance;
    best_metrics.km1 = metrics::km1(_hg);
    best_metrics.cut = metrics::hyperedgeCut(_hg);
    return false;
  }

  void logPQStats() {
    DBG << V(_ordering);
    DBG << V(_inverse_ordering);
    for (const PartitionID& part : _inverse_ordering) {
      DBG << V(part) << V(_pq[PREV].size(part)) << V(_pq[NEXT].size(part)) << V(_hg.partWeight(part));
    }
  }

  void activateMovableHNs() {
    for (const HypernodeID& hn : _hg.nodes()) {
      for (const std::size_t& direction : kDirections) {
        if (_fixtures[direction][hn] == 0 && !_hg.active(hn)) {
          _hg.activate(hn);
        }
      }
    }
  }

  /*********************************************************************************
   * SCHEDULING METHODS
   *********************************************************************************/
  void refreshTopologicalOrdering() {
    DBG << "refreshUnderlyingTopologicalOrdering()";
    _fixtures[PREV].clear();
    _fixtures[NEXT].clear();
    _fixtures[PREV].resize(_hg.initialNumNodes());
    _fixtures[NEXT].resize(_hg.initialNumNodes());
    _pq[PREV].clear();
    _pq[NEXT].clear();
    QuotientGraph<DFSCycleDetector> qg(_hg, _context);
    _ordering.clear();
    const auto& new_order = qg.computeStrictTopologicalOrdering();
    _ordering.insert(_ordering.begin(), new_order.begin(), new_order.end());
    _inverse_ordering.clear();
    _inverse_ordering.resize(_context.partition.k);
    for (PartitionID part = 0; part < _context.partition.k; ++part) {
      _inverse_ordering[_ordering[part]] = part;
    }

    initializeFixtures();
    ASSERT([&]() {
      ASSERT_THAT_FIXTURES_ARE_CONSISTENT();
      return true;
    }());
    processFixtureStateChanges();
  }

  bool orderingStillTopological() const {
    QuotientGraph<DFSCycleDetector> qg(_hg, _context);
    for (PartitionID u = 0; u < _context.partition.k; ++u) {
      for (const auto& edge : qg.outs(u)) {
        const PartitionID& v = edge.first;
        const HyperedgeWeight& weight = edge.second;
        if (weight == 0) { // not a real edge
          continue;
        }
        if (_ordering[u] > _ordering[v]) {
          return false;
        }
      }
    }

    return true;
  }

  std::size_t reverse(const std::size_t direction) const {
    return 1 - direction;
  }

  std::pair<PartitionID, PartitionID> selectPath() const {
    PartitionID start = 0;
    PartitionID end = 0;
    for (PartitionID part = 0; part < _context.partition.k; ++part) {
      if (_hg.partWeight(part) > _hg.partWeight(start)) {
        start = part;
      }
    }
    for (PartitionID part = 0; part < _context.partition.k; ++part) {
      if (_ordering[part] < _ordering[start] && _hg.partWeight(part) <= _hg.partWeight(end)) {
        end = part;
      } else if (_hg.partWeight(part) < _hg.partWeight(end)) {
        end = part;
      }
    }
    ASSERT(start != end, V(start));
    DBG << V(start) << "-->" << V(end);
    return {end, start};

    /*PartitionID last_overloaded_part = Hypergraph::kInvalidPartition;
    PartitionID last_underloaded_part = Hypergraph::kInvalidPartition;
    std::pair<PartitionID, PartitionID> best_path{Hypergraph::kInvalidPartition, Hypergraph::kInvalidPartition};

    HyperedgeWeight best_path_weight = 0;
    HyperedgeWeight current_forward_path_weight = 0; // starts at overloaded part
    HyperedgeWeight current_backward_path_weight = 0; // starts at underloaded part

    for (std::size_t i = 0; i < _inversed_ordering.size(); ++i) {
      const PartitionID part = _inversed_ordering[i];
      if (i > 0) {
        const PartitionID last_part = _inversed_ordering[i - 1];
        if (!_pq[NEXT].empty(last_part)) {
          current_forward_path_weight += _pq[NEXT].maxKey(last_part);
        }
        if (!_pq[PREV].empty(part)) {
          current_backward_path_weight += _pq[PREV].maxKey(part);
        }
      }

      const HyperedgeWeight max_weight = _context.partition.max_part_weights[part];
      const bool is_overloaded = _hg.partSize(part) > (1.0 + _context.partition.epsilon) * max_weight;
      const bool is_underloaded = _hg.partSize(part) < (1.0 - _context.partition.epsilon) * max_weight;

      if (is_overloaded) {
        if (last_underloaded_part != Hypergraph::kInvalidPartition && best_path_weight < current_backward_path_weight) {
          best_path.first = part;
          best_path.second = last_underloaded_part;
          best_path_weight = current_backward_path_weight;
        }
        last_overloaded_part = part;
        current_forward_path_weight = 0;
      }
      if (is_underloaded) {
        if (last_overloaded_part != Hypergraph::kInvalidPartition && best_path_weight < current_forward_path_weight) {
          best_path.first = last_overloaded_part;
          best_path.second = part;
          best_path_weight = current_forward_path_weight;
        }
        last_underloaded_part = part;
        current_backward_path_weight = 0;
      }
    }
    return {6, 1};
    //ASSERT(best_path.first != Hypergraph::kInvalidPartition && best_path.second != Hypergraph::kInvalidPartition);
    //return best_path;*/
  }

  /*********************************************************************************
   * GAIN UPDATE METHODS
   *********************************************************************************/
  void move(const HypernodeID moved_hn, const PartitionID from_part, const PartitionID to_part) {
    DBG << "move(" << V(moved_hn) << V(from_part) << V(to_part) << ")";

    _hg.mark(moved_hn);
    _new_adjacent_part.resetUsedEntries();

    ASSERT(!_gain_cache.entryExists(moved_hn, from_part), V(moved_hn) << V(from_part));

    if (!_gain_cache.entryExists(moved_hn, to_part)) {
      _gain_cache.initializeEntry(moved_hn, to_part, _nonadjacent_gain_cache[moved_hn]);
    }

    bool moved_hn_remains_connected_to_from_part = false;
    for (const HyperedgeID& he : _hg.incidentEdges(moved_hn)) {
      const HypernodeID pins_in_source_part_after = _hg.pinCountInPart(he, from_part);
      moved_hn_remains_connected_to_from_part |= pins_in_source_part_after != 0;

      if (pins_in_source_part_after == 0 && _hg.pinCountInPart(he, to_part) != 1) {
        ASSERT(_gain_cache.entryExists(moved_hn), V(moved_hn));
        for (const PartitionID& part : _gain_cache.adjacentParts(moved_hn)) {
          if (part != from_part && part != to_part) {
            _gain_cache.updateExistingEntry(moved_hn, part, -_hg.edgeWeight(he));
          }
        }
      } else if (pins_in_source_part_after != 0 && _hg.pinCountInPart(he, to_part) == 1) {
        ASSERT(_gain_cache.entryExists(moved_hn), V(moved_hn));
        for (const PartitionID& part : _gain_cache.adjacentParts(moved_hn)) {
          if (part != from_part && part != to_part) {
            _gain_cache.updateExistingEntry(moved_hn, part, _hg.edgeWeight(he));
          }
        }
      }

      updateNeighborsOfMovedHN(moved_hn, from_part, to_part, he);
    }
    _gain_cache.updateFromAndToPartOfMovedHN(moved_hn, from_part, to_part, moved_hn_remains_connected_to_from_part);

    updateFixtures(moved_hn, from_part, to_part);
    processFixtureStateChanges();
  }

  void updateNeighborsOfMovedHN(const HypernodeID moved_hn,
                                const PartitionID from_part,
                                const PartitionID to_part,
                                const HyperedgeID he) {
    const HypernodeID pin_count_from_part_before_move = _hg.pinCountInPart(he, from_part) + 1;
    const HypernodeID pin_count_to_part_after_move = _hg.pinCountInPart(he, to_part);
    const bool move_decreased_connectivity = pin_count_from_part_before_move - 1 == 0;
    const bool move_increased_connectivity = pin_count_to_part_after_move == 1;
    const HyperedgeWeight he_weight = _hg.edgeWeight(he);
    const PinState pin_state(pin_count_from_part_before_move == 1,
                             pin_count_to_part_after_move == 1,
                             pin_count_from_part_before_move == 2,
                             pin_count_to_part_after_move == 2);

    for (const HypernodeID& pin : _hg.pins(he)) {
      DBGC(pin == hn_to_debug) << "updateNeighborsOfMovedHN.for(" << V(pin) << V(_hg.partID(pin))
                               << V(moved_hn) << V(from_part) << V(to_part) << V(he) << ")";
      if (!_hg.marked(pin) && _hg.active(pin)) {
        ASSERT(pin != moved_hn, V(pin));
        connectivityUpdateForCache(pin, from_part, to_part, he,
                                   move_decreased_connectivity,
                                   move_increased_connectivity);
        deltaGainUpdates<true>(pin, from_part, to_part, he);
      } else if (pin != moved_hn) {
        connectivityUpdateForCache(pin, from_part, to_part, he,
                                   move_decreased_connectivity,
                                   move_increased_connectivity);
        deltaGainUpdates<false>(pin, from_part, to_part, he);
      }
    }
  }

  // deltaGainUpdate for pin after another pin in he was moved from from_part to to_part
  template<bool update_pq>
  void deltaGainUpdates(const HypernodeID pin,
                        const PartitionID from_part,
                        const PartitionID to_part,
                        const HyperedgeID he) {
    const HypernodeID pin_count_from_part_before_move = _hg.pinCountInPart(he, from_part) + 1;
    const HypernodeID pin_count_to_part_after_move = _hg.pinCountInPart(he, to_part);
    const HyperedgeWeight he_weight = _hg.edgeWeight(he);
    const PinState pin_state(pin_count_from_part_before_move == 1,
                             pin_count_to_part_after_move == 1,
                             pin_count_from_part_before_move == 2,
                             pin_count_to_part_after_move == 2);
    const PartitionID source_part = _hg.partID(pin);

    DBGC(pin == hn_to_debug)
    << "deltaGainUpdates(" << V(pin) << V(_hg.partID(pin)) << V(from_part) << V(to_part) << V(he) << ")";

    if (source_part == from_part) {
      if (pin_state.two_pins_in_from_part_before) {
        for (const PartitionID& part : _gain_cache.adjacentParts(pin)) {
          if (_new_adjacent_part.get(pin) != part) {
            if (update_pq) {
              updatePin(pin, part, he, he_weight);
            }
            _gain_cache.updateExistingEntry(pin, part, he_weight);
          }
        }
      }
    } else if (source_part == to_part && pin_state.two_pins_in_to_part_after) {
      for (const PartitionID& part : _gain_cache.adjacentParts(pin)) {
        if (_new_adjacent_part.get(pin) != part) {
          if (update_pq) {
            updatePin(pin, part, he, -he_weight);
          }
          _gain_cache.updateExistingEntry(pin, part, -he_weight);
        }
      }
    }

    if (pin_state.one_pin_in_from_part_before && _gain_cache.entryExists(pin, from_part)) {
      if (update_pq) {
        updatePin(pin, from_part, he, -he_weight);
      }
      _gain_cache.updateExistingEntry(pin, from_part, -he_weight);
    }

    if (pin_state.one_pin_in_to_part_after && _new_adjacent_part.get(pin) != to_part) {
      if (update_pq) {
        updatePin(pin, to_part, he, he_weight);
      }
      _gain_cache.updateExistingEntry(pin, to_part, he_weight);
    }
  }

  void updatePin(const HypernodeID pin, const PartitionID part, const HyperedgeID he, const Gain delta) {
    const PartitionID source_part = _hg.partID(pin);
    if (!isAdjacent(source_part, part)) {
      return;
    }
    const std::size_t direction = directionToAdjacent(_hg.partID(pin), part);
    if (_fixtures[direction][pin] > 0) {
      return;
    }
    ASSERT(_pq[direction].contains(pin, source_part));
    _pq[direction].updateKeyBy(pin, source_part, delta);
  }

  std::size_t directionToAdjacent(const PartitionID from, const PartitionID to) {
    ASSERT(isAdjacent(from, to), V(from) << V(to) << V(_ordering));
    return adjacentPart(from, PREV) == to ? PREV : NEXT;
  }

  void connectivityUpdateForCache(const HypernodeID pin, const PartitionID from_part,
                                  const PartitionID to_part, const HyperedgeID he,
                                  const bool move_decreased_connectivity,
                                  const bool move_increased_connectivity) {
    if (move_decreased_connectivity && _gain_cache.entryExists(pin, from_part) &&
        !hypernodeIsConnectedToPart(pin, from_part)) {
      _gain_cache.removeEntryDueToConnectivityDecrease(pin, from_part);
    }
    if (move_increased_connectivity && !_gain_cache.entryExists(pin, to_part)) {
      _gain_cache.addEntryDueToConnectivityIncrease(pin, to_part, gainInducedByHypergraph(pin, to_part));
      _new_adjacent_part.set(pin, to_part);
    }
  }

  void connectivityUpdate(const HypernodeID pin, const PartitionID from_part,
                          const PartitionID to_part, const HyperedgeID he,
                          const bool move_decreased_connectivity,
                          const bool move_increased_connectivity) KAHYPAR_ATTRIBUTE_ALWAYS_INLINE {
    ONLYDEBUG(he);
    const PartitionID source_part = _hg.partID(pin);
    if (move_decreased_connectivity && _gain_cache.entryExists(pin, from_part) &&
        !hypernodeIsConnectedToPart(pin, from_part)) {
      if (isAdjacent(_hg.partID(pin), from_part)) {
        const std::size_t direction = directionToAdjacent(_hg.partID(pin), from_part);
        ASSERT(_pq[direction].contains(pin, source_part));
        _pq[direction].updateKey(pin, _hg.partID(pin), _nonadjacent_gain_cache[pin]);
      }
      _gain_cache.removeEntryDueToConnectivityDecrease(pin, from_part);
    }
    if (move_increased_connectivity && !_gain_cache.entryExists(pin, to_part)) {
      Gain gain = GainCache::kNotCached;
      if (_gain_cache.entryExists(pin, to_part)) {
        gain = _gain_cache.entry(pin, to_part);
      } else {
        gain = gainInducedByHypergraph(pin, to_part);
        _gain_cache.addEntryDueToConnectivityIncrease(pin, to_part, gain);
      }
      if (likely(!_hg.isFixedVertex(pin)) && isAdjacent(_hg.partID(pin), to_part)) {
        const std::size_t direction = directionToAdjacent(_hg.partID(pin), to_part);
        ASSERT(_pq[direction].contains(pin, source_part));
        _pq[direction].updateKey(pin, _hg.partID(pin), gain);
      }
      _new_adjacent_part.set(pin, to_part);
    }
  }

  void enableIfDisabled(const HypernodeID hn) {
    DBG << "enableIfDisabled(" << V(hn) << ")";
    const PartitionID part = _hg.partID(hn);
    const bool inconsistent_in_prev = (_fixtures[PREV][hn] == 0 && !_pq[PREV].contains(hn, part));
    const bool inconsistent_in_next = (_fixtures[NEXT][hn] == 0 && !_pq[NEXT].contains(hn, part));
    if (!inconsistent_in_next && !inconsistent_in_prev) { // must check both for first / last part
      return;
    }
    DBGC(hn == hn_to_debug) << "Enabling disabled" << V(hn);

    for (const HyperedgeID& he : _hg.incidentHeadEdges(hn)) {
      for (const HypernodeID& tail : _hg.tails(he)) {
        if (part == _hg.partID(tail)) {
          ++_fixtures[NEXT][tail];
          ++_fixtures[PREV][hn];
          if (_fixtures[NEXT][tail] == 1) {
            _state_changes[NEXT].set(tail, DEACTIVATE);
          }
        }
      }
    }
    for (const HyperedgeID& he : _hg.incidentTailEdges(hn)) {
      for (const HypernodeID& head : _hg.heads(he)) {
        if (part == _hg.partID(head)) {
          ++_fixtures[PREV][head];
          ++_fixtures[NEXT][hn];
          if (_fixtures[PREV][head] == 1) {
            _state_changes[PREV].set(head, DEACTIVATE);
          }
        }
      }
    }

    if (_fixtures[PREV][hn] == 0) {
      DBGC(hn == hn_to_debug) << "_state_change in PREV to ACTIVE" << V(hn);
      _state_changes[PREV].set(hn, ACTIVATE);
    }
    if (_fixtures[NEXT][hn] == 0) {
      DBGC(hn == hn_to_debug) << "_state_change in NEXT to ACTIVE" << V(hn);
      _state_changes[NEXT].set(hn, ACTIVATE);
    }
  }

  void updateFixtures(const HypernodeID moved_hn, const PartitionID from_part, const PartitionID to_part) {
    DBGC(moved_hn == hn_to_debug || moved_hn == 1803) << "updateFixtures(" << V(moved_hn) << ","
                                                      << V(from_part) << "," << V(to_part) << ")";
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
    return _ordering[part] == 0;
  }

  bool isLastPartition(const PartitionID part) const {
    return _ordering[part] == _context.partition.k - 1;
  }

  bool isAdjacent(const PartitionID a, const PartitionID b) const {
    return std::abs(_ordering[a] - _ordering[b]) == 1;
  }

  bool isValidDirection(const PartitionID part, std::size_t direction) const {
    return !((direction == NEXT && isLastPartition(part)) || (direction == PREV && isFirstPartition(part)));
  }

  Gain cachedGainValue(const HypernodeID hn, const PartitionID to) const {
    DBGC(hn == hn_to_debug) << V(hn) << V(to) << V(_gain_cache.entryExists(hn, to)) << V(_nonadjacent_gain_cache[hn]);
    return _gain_cache.entryExists(hn, to) ? _gain_cache.entry(hn, to)
                                           : _nonadjacent_gain_cache[hn]; // TODO entryExists really necessary?
  }

  void processFixtureStateChanges() {
    DBG << "processFixtureStateChanges(): " << V(_state_changes[PREV].size()) << V(_state_changes[NEXT].size());
    for (const std::size_t& direction : kDirections) {
      for (const auto& hn : _state_changes[direction].usedEntries()) {
        const PartitionID part = _hg.partID(hn);
        if (!isValidDirection(part, direction)) {
          DBGC(hn == hn_to_debug) << "Invalid direction" << V(part) << V(direction);
          continue;
        }

        if (_state_changes[direction].get(hn) == ACTIVATE) {
          const PartitionID to_part = adjacentPart(part, direction);
          DBGC(hn == hn_to_debug) << "Activate HN" << V(hn) << V(direction) << V(part) << V(to_part)
                                  << V(cachedGainValue(hn, to_part));
          if (_pq[direction].contains(hn, part)) {
            _pq[direction].updateKey(hn, part, cachedGainValue(hn, to_part));
          } else {
            _pq[direction].insert(hn, part, cachedGainValue(hn, to_part));
          }
          _pq[direction].enablePart(part);
          if (_hg.marked(hn)) {
            _gain_cache.clear(hn);
            initializeGainCacheFor(hn);
            _hg.unmark(hn);
          }
          if (!_hg.active(hn)) {
            DBGC(hn == hn_to_debug) << "_hg.activate(" << V(hn) << ")";
            _hg.activate(hn);
          }
        } else if (_state_changes[direction].get(hn) == DEACTIVATE) {
          DBGC(hn == hn_to_debug) << "Deactivate HN" << V(hn) << V(direction) << V(part);
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

  void ASSERT_THAT_FIXTURES_ARE_CONSISTENT() const {
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

  Gain gainInducedByHyperedge(const HypernodeID hn, const HyperedgeID he, const PartitionID target_part) const {
    const HypernodeID pins_in_source_part = _hg.pinCountInPart(he, _hg.partID(hn));
    const HypernodeID pins_in_target_part = _hg.pinCountInPart(he, target_part);
    const HyperedgeWeight he_weight = _hg.edgeWeight(he);
    Gain gain = pins_in_source_part == 1 ? he_weight : 0;
    gain -= pins_in_target_part == 0 ? he_weight : 0;
    return gain;
  }

  Gain gainInducedByHypergraph(const HypernodeID hn, const PartitionID target_part) const {
    ASSERT(target_part != _hg.partID(hn), V(hn) << V(target_part));
    Gain gain = 0;
    for (const HyperedgeID& he : _hg.incidentEdges(hn)) {
      ASSERT(_hg.edgeSize(he) > 1, V(he));
      gain += gainInducedByHyperedge(hn, he, target_part);
    }
    return gain;
  }

  PartitionID adjacentPart(const PartitionID part, const std::size_t index) const {
    if (index == PREV) {
      return _ordering[part] == 0 ? Hypergraph::kInvalidPartition : _inverse_ordering[_ordering[part] - 1];
    } else {
      ASSERT(index == NEXT);
      return _ordering[part] == _context.partition.k - 1 ? Hypergraph::kInvalidPartition : _inverse_ordering[
        _ordering[part] + 1];
    }
  }

  void initializeGainCache() {
    for (const HypernodeID& hn : _hg.nodes()) {
      initializeGainCacheFor(hn);
    }
  }

  void initializeGainCacheFor(const HypernodeID hn) {
    DBGC(hn == hn_to_debug) << "initializeGainCacheFor(" << V(hn) << ")";

    _tmp_gains.clear();
    const PartitionID source_part = _hg.partID(hn);
    HyperedgeWeight internal = 0;
    for (const HyperedgeID& he : _hg.incidentEdges(hn)) {
      const HyperedgeWeight he_weight = _hg.edgeWeight(he);
      internal += _hg.pinCountInPart(he, source_part) != 1 ? he_weight : 0;
      for (const PartitionID& part : _hg.connectivitySet(he)) {
        _tmp_gains[part] += he_weight;
      }
    }

    for (const auto& target_part : _tmp_gains) {
      if (target_part.key == source_part) {
        ASSERT(!_gain_cache.entryExists(hn, source_part), V(hn) << V(source_part));
        continue;
      }
      ASSERT(target_part.value - internal == gainInducedByHypergraph(hn, target_part.key),
             V(gainInducedByHypergraph(hn, target_part.key)) << V(target_part.value - internal));
      DBGC(hn == hn_to_debug) << V(target_part.key) << V(target_part.value - internal);
      _gain_cache.initializeEntry(hn, target_part.key, target_part.value - internal);
    }

    _nonadjacent_gain_cache[hn] = -internal;
  }

  bool hypernodeIsConnectedToPart(const HypernodeID pin, const PartitionID part) const {
    for (const HyperedgeID& he : _hg.incidentEdges(pin)) {
      if (_hg.pinCountInPart(he, part) > 0) {
        return true;
      }
    }
    return false;
  }

  void ASSERT_THAT_GAIN_CACHE_IS_VALID() {
    ASSERT([&]() {
      for (const HypernodeID& hn : _hg.nodes()) {
        ASSERT_THAT_CACHE_IS_VALID_FOR_HN(hn);
      }
      return true;
    }(), "Gain Cache inconsistent");
  }

  void ASSERT_THAT_CACHE_IS_VALID_FOR_HN(const HypernodeID hn) const {
    if (_gain_cache.entryExists(hn)) {
      std::vector<bool> adjacent_parts(_context.partition.k, false);
      for (PartitionID part = 0; part < _context.partition.k; ++part) {
        if (hypernodeIsConnectedToPart(hn, part)) {
          adjacent_parts[part] = true;
        }
        if (_gain_cache.entry(hn, part) != GainCache::kNotCached) {
          ASSERT(_gain_cache.entryExists(hn, part), V(hn) << V(part));
          ASSERT(_gain_cache.entry(hn, part) == gainInducedByHypergraph(hn, part),
                 V(hn) << V(part) << V(_gain_cache.entry(hn, part)) <<
                       V(gainInducedByHypergraph(hn, part)) << V(_hg.partID(hn)));
          ASSERT(hypernodeIsConnectedToPart(hn, part), V(hn) << V(part));
        } else if (_hg.partID(hn) != part && !hypernodeIsConnectedToPart(hn, part)) {
          ASSERT(!_gain_cache.entryExists(hn, part), V(hn) << V(part)
                                                           << "_hg.partID(hn) != part");
          ASSERT(_gain_cache.entry(hn, part) == GainCache::kNotCached, V(hn) << V(part));
        }
        if (_hg.partID(hn) == part) {
          ASSERT(!_gain_cache.entryExists(hn, part), V(hn) << V(part)
                                                           << "_hg.partID(hn) == part");
          ASSERT(_gain_cache.entry(hn, part) == GainCache::kNotCached,
                 V(hn) << V(part));
        }
      }
      for (const PartitionID& part : _gain_cache.adjacentParts(hn)) {
        ASSERT(adjacent_parts[part], V(part));
      }
    }
  }

  Hypergraph& _hg;
  Context _context;
  std::array<KWayRefinementPQ, TARGETS> _pq;
  std::array<std::vector<HypernodeID>, TARGETS> _fixtures;
  std::array<ds::FastResetArray<std::size_t>, 2> _state_changes;
  std::vector<PartitionID> _ordering;
  std::vector<PartitionID> _inverse_ordering;
  GainCache _gain_cache;
  ds::SparseMap<PartitionID, Gain> _tmp_gains;
  std::vector<Gain> _nonadjacent_gain_cache;
  ds::FastResetArray<PartitionID> _new_adjacent_part;
  std::vector<Move> _moves;
  std::vector<HypernodeID> _moved_hns;
};

constexpr std::array<std::size_t, 2> AcyclicHardRebalanceRefiner::kDirections;
} // namespace kahypar