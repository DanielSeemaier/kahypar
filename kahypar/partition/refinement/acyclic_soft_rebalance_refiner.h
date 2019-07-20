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

  AcyclicSoftRebalanceRefiner(Hypergraph& hypergraph, const Context& context, AdjacencyMatrixQuotientGraph<DFSCycleDetector>& qg) :
    _hg(hypergraph),
    _context(context),
    _qg(qg),
    _pq(context.partition.k),
    _min_successor_part(_hg.initialNumNodes()),
    _max_predecessor_part(_hg.initialNumNodes()),
    _updated_neighbors(_hg.initialNumNodes(), false) {}

  ~AcyclicSoftRebalanceRefiner() override = default;

  AcyclicSoftRebalanceRefiner(const AcyclicSoftRebalanceRefiner&) = delete;
  AcyclicSoftRebalanceRefiner& operator=(const AcyclicSoftRebalanceRefiner&) = delete;

  AcyclicSoftRebalanceRefiner(AcyclicSoftRebalanceRefiner&&) = delete;
  AcyclicSoftRebalanceRefiner& operator=(AcyclicSoftRebalanceRefiner&&) = delete;

  const std::vector<Move>& moves() const {
    return _moves;
  }

  void printSummary() const override {

  }

  void performMovesAfterHardRebalance(const std::vector<Move>& moves,
                                      std::vector<HypernodeID>& refinement_nodes,
                                      const UncontractionGainChanges& changes) {
    ASSERT(orderingStillTopological(), V(orderingStillTopologicalDebug()));

    _hg.resetHypernodeState();
    for (std::size_t i = moves.size(); i > 0; --i) {
      const Move& move = moves[i - 1];
      _hg.changeNodePart(move.hn, move.to, move.from);
    }

    for (const Move& move : moves) {
      if (!_hg.active(move.hn)) {
        _hg.activate(move.hn);
      }
      _hg.mark(move.hn);

      _hg.changeNodePart(move.hn, move.from, move.to);
      if (!_skip_pq_updates && !_pq_inst[move.hn].empty()) {
        ASSERT(_pq.contains(move.hn, move.from));
        _pq.remove(move.hn, move.from);
        _pq_inst[move.hn].clear();
      }
      this->move(move.hn, move.from, move.to);

      _hg.unmark(move.hn);
      resetMinMaxNeighborsFor(move.hn);
      initializeMinMaxNeighborsFor(move.hn);
    }

    _gain_cache.resetDelta();
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
      const bool success = _qg.testAndUpdateBeforeMovement(move.hn, move.to, move.from, changed);
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
      const bool success = _qg.testAndUpdateBeforeMovement(move.hn, move.from, move.to, changed);
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
      const bool success = _qg.testAndUpdateBeforeMovement(hn, from_part, to_part);
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
    ASSERT(AdjacencyMatrixQuotientGraph<DFSCycleDetector>(_hg, _context).isAcyclic(),
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
  }

  /*********************************************************************************
   * SCHEDULING METHODS
   *********************************************************************************/
  void refreshTopologicalOrdering() {
    DBG << "refreshUnderlyingTopologicalOrdering()";
    _skip_pq_updates = false;

    _hg.resetHypernodeState();
    _ordering.clear();
    const auto& new_order = _qg.topologicalOrdering();
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

  /*********************************************************************************
   * GAIN UPDATE METHODS
   *********************************************************************************/
  void move(const HypernodeID moved_hn, const PartitionID from_part, const PartitionID to_part) {
    DBG << "move(" << V(moved_hn) << V(from_part) << V(to_part) << ")";
    ASSERT(_hg.marked(moved_hn));

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


  void updatePin(const HypernodeID pin, const PartitionID part, const HyperedgeID he, const Gain delta) {
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

#ifdef KAHYPAR_USE_ASSERTIONS
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

  bool isConnectedTo(const HypernodeID pin, const PartitionID part) const {
    for (const HyperedgeID& he : _hg.incidentEdges(pin)) {
      if (_hg.pinCountInPart(he, part) > 0) {
        return true;
      }
    }
    return false;
  }
#endif // KAHYPAR_USE_ASSERTIONS

  Hypergraph& _hg;
  const Context& _context;
  AdjacencyMatrixQuotientGraph<DFSCycleDetector>& _qg;

  KWayRefinementPQ _pq;
  std::vector<BinaryMaxHeap<PartitionID, Gain>> _pq_inst;

  std::vector<PartitionID> _min_successor_part;
  std::vector<PartitionID> _max_predecessor_part;

  std::vector<Move> _moves;

  ds::FastResetArray<bool> _updated_neighbors;

  std::size_t _num_moves = 0;
  std::size_t _num_touched_hns = 0;
  std::size_t _num_rejected_due_to_balance = 0;
  std::size_t _num_rejected_due_to_km1 = 0;

  double _improved_imbalance_by = 0.0;
};
} // namespace kahypar