#pragma once

namespace kahypar {
class AcyclicSoftRebalanceRefiner final : public IRefiner {
 private:
  static constexpr bool debug = false;
  static constexpr HypernodeID hn_to_debug = 222;

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
  using GainCache = KwayGainCache<Gain>;
  using KWayRefinementPQ = ds::KWayPriorityQueue<HypernodeID, Gain, std::numeric_limits<Gain>>;

  AcyclicSoftRebalanceRefiner(Hypergraph& hypergraph, const Context& context, QuotientGraph<DFSCycleDetector>& qg) :
    _hg(hypergraph),
    _context(context),
    _qg(qg),
    _pq(context.partition.k),
    _new_adjacent_part(_hg.initialNumNodes(), Hypergraph::kInvalidPartition),
    _min_successor_part(_hg.initialNumNodes()),
    _max_predecessor_part(_hg.initialNumNodes()),
    _gain_cache(_hg.initialNumNodes(), _context.partition.k),
    _tmp_gains(_context.partition.k, 0),
    _nonadjacent_gain_cache(_hg.initialNumNodes()),
    _updated_neighbors(_hg.initialNumNodes(), false) {}

  ~AcyclicSoftRebalanceRefiner() override = default;

  AcyclicSoftRebalanceRefiner(const AcyclicSoftRebalanceRefiner&) = delete;
  AcyclicSoftRebalanceRefiner& operator=(const AcyclicSoftRebalanceRefiner&) = delete;

  AcyclicSoftRebalanceRefiner(AcyclicSoftRebalanceRefiner&&) = delete;
  AcyclicSoftRebalanceRefiner& operator=(AcyclicSoftRebalanceRefiner&&) = delete;

  void changeEpsilon(const double epsilon) {
    _context.partition.epsilon = epsilon;
    _context.setupPartWeights(_hg.totalWeight());
  }

  const std::vector<Move>& moves() {
    return _moves;
  }

  std::size_t numMoves() const {
    return _num_moves;
  }

  double improvedImbalanceBy() const {
    return _improved_imbalance_by;
  }

  std::size_t numRejectedDueToBalance() const {
    return _num_rejected_due_to_balance;
  }

  std::size_t numRejectedDueToKM1() const {
    return _num_rejected_due_to_km1;
  }

  void enableRefinementNodes(std::vector<HypernodeID>& refinement_nodes) {
    LOG << "enableRefinementNodes(" << refinement_nodes << ")";
    _hg.resetHypernodeState();
    _updated_neighbors.resetUsedEntries();

    for (const HypernodeID& hn : refinement_nodes) {
      DBGC(hn == hn_to_debug) << "Reset" << V(hn);
      _gain_cache.clear(hn);
      initializeGainCacheFor(hn);
      reinitNodeWithNeighbors(hn, false);
    }

    _gain_cache.resetDelta();
  }

  std::size_t numPosGains() const {
    return _num_pos_gain_hns;
  }

  std::size_t numNegGains() const {
    return _num_neg_gain_hns;
  }

  std::size_t numZeroGains() const {
    return _num_zero_gain_hns;
  }

 private:
  void initializeImpl(const HyperedgeWeight max_gain) final {
    LOG << "initializeImpl()";
    if (!_is_initialized) {
      _pq.initialize(_hg.initialNumNodes());
      for (HypernodeID hn = 0; hn < _hg.initialNumNodes(); ++hn) {
        _pq_inst.emplace_back(_context.partition.k);
      }
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
    _hg.resetHypernodeState();
    if (!_skip_pq_updates) {
      _skip_pq_updates = !orderingStillTopological();
    }

    for (std::size_t i = moves.size(); i > 0; --i) {
      const Move& move = moves[i - 1];
      bool changed;
      const bool success = _qg.update(move.hn, move.to, move.from, changed);
      ASSERT(success);
      if (changed) {
        _skip_pq_updates = true;
      }
      _hg.changeNodePart(move.hn, move.to, move.from);
    }

    if (!_skip_pq_updates) {
      _skip_pq_updates = !orderingStillTopological();
    }

    for (const Move& move : moves) {
      DBG << "performMovesAndUpdateCacheImpl.move(" << V(move.hn) << V(move.from) << V(move.to) << ")";
      if (!_hg.active(move.hn)) {
        _hg.activate(move.hn);
      }
      _hg.mark(move.hn);
      bool changed = false;
      const bool success = _qg.update(move.hn, move.from, move.to, changed);
      ASSERT(success);
      if (changed) {
        _skip_pq_updates = true;
      }
      _hg.changeNodePart(move.hn, move.from, move.to);
      if (!_skip_pq_updates && !_pq_inst[move.hn].empty()) {
        ASSERT(_pq.contains(move.hn, move.from));
        _pq.remove(move.hn, move.from);
        _pq_inst[move.hn].clear();
      }
      this->move(move.hn, move.from, move.to);
    }

    _gain_cache.resetDelta();
  }

  bool refineImpl(std::vector<HypernodeID>& refinement_nodes,
                  const std::array<HypernodeWeight, 2>&,
                  const UncontractionGainChanges&,
                  Metrics& best_metrics) final {
    DBG << "refineImpl(" << refinement_nodes << "): imbalance" << best_metrics.imbalance;
    ASSERT(best_metrics.imbalance == metrics::imbalance(_hg, _context));
    _moves.clear();
    _hg.resetHypernodeState();

    if (!orderingStillTopological() || _skip_pq_updates) {
      refreshTopologicalOrdering();
      LOG << "refreshed";
    }
    ASSERT(orderingStillTopological());
    ASSERT(orderingStillTopologicalDebug());

    ASSERT([&]() {
      for (const HypernodeID& hn : _hg.nodes()) {
        if (_hg.marked(hn)) {
          continue;
        }

        if (_max_predecessor_part[hn] != _min_successor_part[hn]) {
          ASSERT(_pq.contains(hn, _hg.partID(hn)));
          ASSERT(!_pq_inst[hn].empty());
        } else {
          ASSERT(!_pq.contains(hn, _hg.partID(hn)));
          ASSERT(_pq_inst[hn].empty());
        }
      }
      return true;
    }());

    std::size_t touched_hns_since_improvement = 0;
    int last_accepted_index = -1;
    const HyperedgeWeight initial_km1 = best_metrics.km1;
    const double initial_imbalance = best_metrics.imbalance;
    HyperedgeWeight current_km1 = best_metrics.km1;

    for (PartitionID part = 0; part < _context.partition.k; ++part) {
      if (isOverloaded(part)) {
        _pq.enablePart(part);
      } else {
        _pq.disablePart(part);
      }
    }

    std::vector<HypernodeID> marked_hns;

    while (touched_hns_since_improvement < 25 && !_pq.empty()) { // TODO make that a config parameter
      HypernodeID max_gain_hn;
      PartitionID from_part;
      Gain gain;
      _pq.deleteMax(max_gain_hn, gain, from_part);
      if (gain > 0) {
        ++_num_pos_gain_hns;
      } else if (gain < 0) {
        ++_num_neg_gain_hns;
      } else {
        ++_num_zero_gain_hns;
      }
      ASSERT(from_part == _hg.partID(max_gain_hn));
      ASSERT(!_hg.marked(max_gain_hn));
      const PartitionID to_part = _pq_inst[max_gain_hn].top();
      ASSERT(_max_predecessor_part[max_gain_hn] <= _ordering[to_part]);
      ASSERT(_min_successor_part[max_gain_hn] >= _ordering[from_part]);
      DBG << "Removed" << V(max_gain_hn) << V(gain) << V(from_part) << V(to_part) << V(_hg.partID(max_gain_hn));
      _pq_inst[max_gain_hn].clear();

      if (!_hg.active(max_gain_hn)) {
        _hg.activate(max_gain_hn);
      }
      _hg.mark(max_gain_hn);
      marked_hns.push_back(max_gain_hn);

      ++_num_touched_hns;
      ++touched_hns_since_improvement;

      if (_hg.partWeight(to_part) + _hg.nodeWeight(max_gain_hn) < _hg.partWeight(from_part)) {
        HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
        const bool success = _qg.update(max_gain_hn, from_part, to_part);
        ASSERT(success);
        HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
        Timer::instance().add(_context, Timepoint::cycle_detector,
                              std::chrono::duration<double>(end - start).count());
        _hg.changeNodePart(max_gain_hn, from_part, to_part);
        move(max_gain_hn, from_part, to_part);

        if (!isOverloaded(from_part)) {
          _pq.disablePart(from_part);
        }
        if (!_pq.isEnabled(to_part) && isOverloaded(to_part)) {
          _pq.enablePart(to_part);
        }

        current_km1 -= gain;
        _moves.emplace_back(max_gain_hn, from_part, to_part);
        ++_num_moves;

        if (current_km1 <= initial_km1) {
          touched_hns_since_improvement = 0;
          last_accepted_index = static_cast<int>(_moves.size()) - 1;
          _gain_cache.resetDelta();
        }
      } else {
        ++_num_rejected_due_to_balance;
      }
    }

    ASSERT_THAT_GAIN_CACHE_IS_VALID();

    _gain_cache.rollbackDelta();

    int last_index = static_cast<int>(_moves.size()) - 1;
    DBG << V(last_accepted_index) << V(last_index);
    while (last_index != last_accepted_index) {
      const HypernodeID hn = _moves[last_index].hn;
      const PartitionID from_part = _moves[last_index].to;
      const PartitionID to_part = _moves[last_index].from;
      DBG << "Revert" << V(hn) << V(from_part) << V(to_part);
      _hg.changeNodePart(hn, from_part, to_part);
      HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
      const bool success = _qg.update(hn, from_part, to_part);
      ASSERT(success);
      HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
      Timer::instance().add(_context, Timepoint::cycle_detector,
                            std::chrono::duration<double>(end - start).count());
      _moves.pop_back();
      --last_index;
      --_num_moves;
      ++_num_rejected_due_to_km1;

      // rollbackDelta() does not account for clears, but we currently clear and re-init a cache entry after the node was
      // moved, thus we also have to do that during rollback ...
      _gain_cache.clear(hn);
      initializeGainCacheFor(hn);


      // moving the node back might change the move targets of neighboring nodes
      reinitNodeWithNeighbors(hn, true); // TODO is there a way around this?
    }

    for (const HypernodeID& marked_hn : marked_hns) {
      _hg.unmark(marked_hn);
      resetMinMaxNeighborsFor(marked_hn);
      initializeMinMaxNeighborsFor(marked_hn);
    }

    best_metrics.imbalance = metrics::imbalance(_hg, _context);

    ASSERT_THAT_GAIN_CACHE_IS_VALID();
    ASSERT(_qg.isAcyclic(), "Rollback produced a cyclic quotient graph!");
    ASSERT(QuotientGraph<DFSCycleDetector>(_hg, _context).isAcyclic(),
           "Rollback produced a cyclic quotient graph not detected by the QG!");

    DBG << "Imbalance before _soft_rebalance::refineImpl():" << initial_imbalance;
    DBG << "Imbalance after _soft_rebalance::refineImpl():" << best_metrics.imbalance << V(_num_moves);
    best_metrics.km1 = metrics::km1(_hg);
    best_metrics.cut = metrics::hyperedgeCut(_hg);
    _improved_imbalance_by += initial_imbalance - best_metrics.imbalance;
    DBG << V(_num_moves) << V(_num_touched_hns) << V(initial_km1) << V(best_metrics.km1);
    return false;
  }

  void reinitNodeWithNeighbors(const HypernodeID hn, bool reset = true) {
    if (_skip_pq_updates) {
      return;
    }
    if (reset) {
      _updated_neighbors.resetUsedEntries();
    }

    initializeMinMaxNeighborsFor(hn);
    _updated_neighbors.set(hn, true);

    for (const HyperedgeID& he : _hg.incidentEdges(hn)) {
      for (const HypernodeID& pin : _hg.pins(he)) {
        if (_updated_neighbors.get(pin)) {
          continue;
        }
        _updated_neighbors.set(pin, true);
        initializeMinMaxNeighborsFor(pin);
      }
    }
  }

  bool isOverloaded(const PartitionID part) const {
    return _hg.partWeight(part) > _context.partition.max_part_weights[part];
  }

  void initializeMinMaxNeighbors() {
    for (const HypernodeID& hn : _hg.nodes()) {
      initializeMinMaxNeighborsFor(hn);
    }
  }

  void resetMinMaxNeighborsFor(const HypernodeID hn) {
    _pq_inst[hn].clear();
    _min_successor_part[hn] = _hg.partID(hn);
    _max_predecessor_part[hn] = _hg.partID(hn);
  }

  void initializeMinMaxNeighborsFor(const HypernodeID hn) {
    if (_hg.marked(hn)) {
      return;
    }

    const PartitionID old_min_successor = _min_successor_part[hn];
    const PartitionID old_max_predecessor = _max_predecessor_part[hn];

    _min_successor_part[hn] = _context.partition.k - 1;
    _max_predecessor_part[hn] = 0;

    for (const HyperedgeID& he : _hg.incidentTailEdges(hn)) {
      for (const HypernodeID& head : _hg.heads(he)) {
        _min_successor_part[hn] = std::min(_ordering[_hg.partID(head)], _min_successor_part[hn]);
      }
    }
    for (const HyperedgeID& he : _hg.incidentHeadEdges(hn)) {
      for (const HypernodeID& tail : _hg.tails(he)) {
        ASSERT(orderingStillTopological());
        ASSERT(orderingStillTopologicalDebug());
        ASSERT(_ordering[_hg.partID(tail)] <= _ordering[_hg.partID(hn)],
               V(_ordering[_hg.partID(tail)]) << V(tail) << V(_ordering[_hg.partID(hn)]) << V(hn));
        _max_predecessor_part[hn] = std::max(_ordering[_hg.partID(tail)], _max_predecessor_part[hn]);
      }
    }

    ASSERT(_max_predecessor_part[hn] <= _ordering[_hg.partID(hn)]);
    ASSERT(_min_successor_part[hn] >= _ordering[_hg.partID(hn)]);

    if (old_min_successor != old_max_predecessor) {
      ASSERT(_pq.contains(hn, _hg.partID(hn)), V(hn) << V(_hg.partID(hn)) << V(_pq_inst[hn].size()));
      ASSERT(!_pq_inst[hn].empty(), V(hn));
      _pq.remove(hn, _hg.partID(hn));
      _pq_inst[hn].clear();
    }
    if (_min_successor_part[hn] != _max_predecessor_part[hn]) {
      for (PartitionID topo_part = _max_predecessor_part[hn]; topo_part <= _min_successor_part[hn]; ++topo_part) {
        const PartitionID part = _inverse_ordering[topo_part];
        if (part == _hg.partID(hn)) {
          continue;
        }
        _pq_inst[hn].push(part, cachedGainValue(hn, part));
      }
      _pq.insert(hn, _hg.partID(hn), _pq_inst[hn].topKey());
    }
  }

  void updateMinMaxNeighborsAfterMove(const HypernodeID moved_hn, const PartitionID from_part,
                                      const PartitionID to_part) {
    _updated_neighbors.resetUsedEntries();
    // otherwise initializeMinMaxNeigborsFor() will try to remove moved_hn from _pq, but the moved hn is already
    // removed from the PQ by the caller of the move() method
    _min_successor_part[moved_hn] = from_part;
    _max_predecessor_part[moved_hn] = from_part;
    initializeMinMaxNeighborsFor(moved_hn);
    _updated_neighbors.set(moved_hn, true);


    for (const HyperedgeID& he : _hg.incidentEdges(moved_hn)) {
      for (const HypernodeID& pin : _hg.pins(he)) {
        if (_updated_neighbors.get(pin) || _hg.marked(pin)) {
          continue;
        }
        _updated_neighbors.set(pin, true);
        initializeMinMaxNeighborsFor(pin);
      }
    }
//    for (const HyperedgeID& he : _hg.incidentTailEdges(moved_hn)) {
//      for (const HypernodeID& head : _hg.heads(he)) {
//        if (_updated_neighbors.get(head) || _hg.marked(head)) {
//          continue;
//        }
//        const bool reinitialize = _max_predecessor_part[head] == _ordering[from_part]
//                                  || _max_predecessor_part[head] < _ordering[to_part];
//        DBGC(head == hn_to_debug) << V(head) << V(reinitialize);
//        if (reinitialize) {
//          initializeMinMaxNeighborsFor(head);
//          _updated_neighbors.set(head, true);
//        }
//      }
//    }
//
//    for (const HyperedgeID& he : _hg.incidentHeadEdges(moved_hn)) {
//      for (const HypernodeID& tail : _hg.tails(he)) {
//        if (_updated_neighbors.get(tail) || _hg.marked(tail)) {
//          continue;
//        }
//        const bool reinitialize = _min_successor_part[tail] == _ordering[from_part]
//                                  || _min_successor_part[tail] > _ordering[to_part];
//        DBGC(tail == hn_to_debug) << V(tail) << V(reinitialize);
//        if (reinitialize) {
//          initializeMinMaxNeighborsFor(tail);
//          _updated_neighbors.set(tail, true);
//        }
//      }
//    }
  }

  /*********************************************************************************
   * SCHEDULING METHODS
   *********************************************************************************/
  void refreshTopologicalOrdering() {
    DBG << "refreshUnderlyingTopologicalOrdering()";
    _skip_pq_updates = false;

    _hg.resetHypernodeState();
    _ordering.clear();
    const auto& new_order = _qg.computeStrictTopologicalOrdering();
    _ordering.insert(_ordering.begin(), new_order.begin(), new_order.end());
    _inverse_ordering.clear();
    _inverse_ordering.resize(_context.partition.k);
    for (PartitionID part = 0; part < _context.partition.k; ++part) {
      _inverse_ordering[_ordering[part]] = part;
    }
    _pq.clear();
    for (const HypernodeID& hn : _hg.nodes()) {
      resetMinMaxNeighborsFor(hn);
    }
    initializeMinMaxNeighbors();
  }

  bool orderingStillTopological() const {
    for (PartitionID u = 0; u < _context.partition.k; ++u) {
      for (const auto& edge : _qg.outs(u)) {
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

  bool orderingStillTopologicalDebug() const {
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

  /*********************************************************************************
   * GAIN UPDATE METHODS
   *********************************************************************************/
  void move(const HypernodeID moved_hn, const PartitionID from_part, const PartitionID to_part) {
    DBG << "move(" << V(moved_hn) << V(from_part) << V(to_part) << ")";
    ASSERT(_hg.marked(moved_hn));

    _new_adjacent_part.resetUsedEntries();

    ASSERT(!_gain_cache.entryExists(moved_hn, from_part), V(moved_hn) << V(from_part));
    if (!_gain_cache.entryExists(moved_hn, to_part)) {
      DBGC(moved_hn == hn_to_debug) << V(to_part) << V(_nonadjacent_gain_cache[moved_hn]);
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

    // TODO why doesn't it work without this?!
    _gain_cache.clear(moved_hn);
    initializeGainCacheFor(moved_hn);

    if (!_skip_pq_updates) {
      updateMinMaxNeighborsAfterMove(moved_hn, from_part, to_part);
      ASSERT([&]() {
        for (const HypernodeID& hn : _hg.nodes()) {
          if (_hg.marked(hn)) {
            continue;
          }

          if (_max_predecessor_part[hn] != _min_successor_part[hn]) {
            ASSERT(_pq.contains(hn, _hg.partID(hn)), V(hn) << V(_hg.partID(hn)));
            ASSERT(!_pq_inst[hn].empty());
          } else {
            ASSERT(!_pq.contains(hn, _hg.partID(hn)));
            ASSERT(_pq_inst[hn].empty());
          }
        }
        return true;
      }());
    }

    ASSERT_THAT_GAIN_CACHE_IS_VALID();
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
      if (!_hg.marked(pin) && !_pq_inst[pin].empty()) {
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
    if (_skip_pq_updates) {
      return;
    }

    const PartitionID source_part = _hg.partID(pin);
    if (movableTo(pin, part)) {
      ASSERT(_pq_inst[pin].contains(part));
      const bool is_top = _pq_inst[pin].top() == part;
      _pq_inst[pin].updateKeyBy(part, delta);
      if (is_top) {
        ASSERT(_pq.contains(pin, source_part));
        _pq.updateKey(pin, source_part, _pq_inst[pin].topKey());
      }
    }
  }

  void updatePinTo(const HypernodeID pin, const PartitionID part, const Gain gain) {
    const PartitionID source_part = _hg.partID(pin);
    if (movableTo(pin, part)) {
      ASSERT(_pq_inst[pin].contains(part));
      const bool is_top = _pq_inst[pin].top() == part;
      _pq_inst[pin].updateKey(part, gain);
      if (is_top) {
        ASSERT(_pq.contains(pin, source_part));
        _pq.updateKey(pin, source_part, _pq_inst[pin].topKey());
      }
    }
  }

  bool movableTo(const HypernodeID hn, const PartitionID part) const {
    return _hg.partID(hn) != part
           && _max_predecessor_part[hn] <= _ordering[part]
           && _ordering[part] <= _min_successor_part[hn];
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
      // switch to _nonadjacent_gain_cache
      if (!_skip_pq_updates && movableTo(pin, from_part)) {
        updatePinTo(pin, from_part, _nonadjacent_gain_cache[pin]);
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
      // switch to _gain_cache
      if (!_skip_pq_updates && likely(!_hg.isFixedVertex(pin)) && movableTo(pin, to_part)) {
        updatePinTo(pin, to_part, gain);
      }
      _new_adjacent_part.set(pin, to_part);
    }
  }

  Gain cachedGainValue(const HypernodeID hn, const PartitionID to) const {
    DBGC(hn == hn_to_debug) << V(hn) << V(to) << V(_gain_cache.entryExists(hn, to)) << V(_nonadjacent_gain_cache[hn]);
    return _gain_cache.entryExists(hn, to) ? _gain_cache.entry(hn, to)
                                           : _nonadjacent_gain_cache[hn]; // TODO entryExists really necessary?
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
  QuotientGraph<DFSCycleDetector>& _qg;

  KWayRefinementPQ _pq;
  std::vector<BinaryMaxHeap<PartitionID, Gain>> _pq_inst;

  std::vector<PartitionID> _min_successor_part;
  std::vector<PartitionID> _max_predecessor_part;

  std::vector<PartitionID> _ordering;
  std::vector<PartitionID> _inverse_ordering;
  GainCache _gain_cache;
  ds::SparseMap<PartitionID, Gain> _tmp_gains;
  std::vector<Gain> _nonadjacent_gain_cache;
  ds::FastResetArray<PartitionID> _new_adjacent_part;
  std::vector<Move> _moves;

  ds::FastResetArray<bool> _updated_neighbors;

  std::size_t _num_moves = 0;
  std::size_t _num_touched_hns = 0;
  std::size_t _num_rejected_due_to_balance = 0;
  std::size_t _num_rejected_due_to_km1 = 0;

  std::vector<HypernodeID> _marked_hns;

  bool _skip_pq_updates = false;
  double _improved_imbalance_by = 0.0;

  std::size_t _num_zero_gain_hns = 0;
  std::size_t _num_pos_gain_hns = 0;
  std::size_t _num_neg_gain_hns = 0;
};
} // namespace kahypar