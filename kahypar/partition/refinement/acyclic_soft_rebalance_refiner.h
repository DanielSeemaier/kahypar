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

 public:
  using KWayRefinementPQ = ds::KWayPriorityQueue<HypernodeID, Gain, std::numeric_limits<Gain>>;

  AcyclicSoftRebalanceRefiner(Hypergraph& hypergraph,
                              const Context& context,
                              AdjacencyMatrixQuotientGraph<DFSCycleDetector>& qg,
                              KMinusOneGainManager& gain_manager) :
    _hg(hypergraph),
    _context(context),
    _qg(qg),
    _gain_manager(gain_manager),
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
    LOG << "[SoftRebalance] Iterations with refresh:" << _num_refreshes;
    LOG << "[SoftRebalance] Iterations without refresh:" << _num_no_refreshes;
    LOG << "[SoftRebalance] Number of moves:" << _num_moves;
    LOG << "[SoftRebalance] Number of moves in the last iteration:" << _num_moves_in_last_iteration;
    LOG << "[SoftRebalance] Positive gain moves:" << _num_positive_gain_moves;
    LOG << "[SoftRebalance] Zero gain moves:" << _num_zero_gain_moves;
    LOG << "[SoftRebalance] Negative gain moves:" << _num_negative_gain_moves;
    LOG << "[SoftRebalance] Improved imbalance by:" << _improved_imbalance;
    LOG << "[SoftRebalance] Improved KM1 by:" << _improved_km1;
  }

  void preUncontraction(const HypernodeID representant) override {
    if (_qg_changed) {
      return;
    }
  }

  void postUncontraction(const HypernodeID representant, const std::vector<HypernodeID>&& partners) override {
    if (_qg_changed) {
      return;
    }
  }

 private:
  void initializeImpl(const HyperedgeWeight max_gain) final {
    if (!_is_initialized) {
      _pq.initialize(_hg.initialNumNodes());
      for (HypernodeID hn = 0; hn < _hg.initialNumNodes(); ++hn) {
        _pq_inst.emplace_back(_context.partition.k);
      }
      _is_initialized = true;
    }
    _hg.resetHypernodeState();
    refreshTopologicalOrdering();

    _num_refreshes = 0;
    _num_no_refreshes = 0;
    _num_moves = 0;
    _num_moves_in_last_iteration = 0;
    _num_zero_gain_moves = 0;
    _num_positive_gain_moves = 0;
    _num_negative_gain_moves = 0;
    _improved_km1 = 0;
    _improved_imbalance = 0.0;
  }

  void performMovesAndUpdateCacheImpl(const std::vector<Move>& moves,
                                      std::vector<HypernodeID>& refinement_nodes,
                                      const UncontractionGainChanges& changes) final {
    _qg.topologicalOrdering();
    if (_qg.changed() || _qg_changed) { // cached topological ordering changed, must reset PQs anyways
      _qg_changed = true;
      return;
    }

    _hg.resetHypernodeState();

    for (std::size_t i = moves.size(); i > 0; --i) {
      const Move& move = moves[i - 1];
      _hg.changeNodePart(move.hn, move.to, move.from);
    }

    for (const Move& move : moves) {
      _hg.changeNodePart(move.hn, move.from, move.to);
      this->move(move.hn, move.from, move.to);
    }
  }

  bool refineImpl(std::vector<HypernodeID>& refinement_nodes,
                  const std::array<HypernodeWeight, 2>&,
                  const UncontractionGainChanges&,
                  Metrics& best_metrics) final {
    ASSERT(best_metrics.imbalance == metrics::imbalance(_hg, _context));
    ASSERT(best_metrics.km1 == metrics::km1(_hg));
    _hg.resetHypernodeState();
    _moves.clear();

    if (_qg_changed) {
      refreshTopologicalOrdering();
      ++_num_refreshes;
    } else {
      ++_num_no_refreshes;
    }

    ASSERT([&]() {
      ASSERT_THAT_PQS_CONTAIN_CORRECT_HYPERNODES();
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
      HypernodeID max_gain_hn = Hypergraph::kInvalidHypernodeID;
      PartitionID from_part = Hypergraph::kInvalidPartition;
      Gain gain;
      _pq.deleteMax(max_gain_hn, gain, from_part);
      if (gain > 0) {
        ++_num_positive_gain_moves;
      } else if (gain < 0) {
        ++_num_negative_gain_moves;
      } else {
        ++_num_zero_gain_moves;
      }
      ASSERT(max_gain_hn != Hypergraph::kInvalidHypernodeID);
      ASSERT(from_part == _hg.partID(max_gain_hn));
      ASSERT(_hg.active(max_gain_hn));
      ASSERT(!_hg.marked(max_gain_hn));

      const auto& ordering = _qg.topologicalOrdering();
      const PartitionID to_part = _pq_inst[max_gain_hn].top();
      ASSERT(_max_predecessor_part[max_gain_hn] <= ordering[to_part]);
      ASSERT(_min_successor_part[max_gain_hn] >= ordering[from_part]);

      _pq_inst[max_gain_hn].clear();
      _hg.mark(max_gain_hn);
      marked_hns.push_back(max_gain_hn);

      ++touched_hns_since_improvement;
      if (_hg.partWeight(to_part) + _hg.nodeWeight(max_gain_hn) < _hg.partWeight(from_part)) {
        ++_num_moves;

        // move hypernode
        HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
        const bool success = _qg.testAndUpdateBeforeMovement(max_gain_hn, from_part, to_part);
        ASSERT(success);
        HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
        Timer::instance().add(_context, Timepoint::cycle_detector,
                              std::chrono::duration<double>(end - start).count());

        _hg.changeNodePart(max_gain_hn, from_part, to_part);
        _gain_manager.updateAfterMovement(max_gain_hn, from_part, to_part);
        move(max_gain_hn, from_part, to_part);

        if (!isOverloaded(from_part)) {
          _pq.disablePart(from_part);
        }
        if (!_pq.isEnabled(to_part) && isOverloaded(to_part)) {
          _pq.enablePart(to_part);
        }

        current_km1 -= gain;
        ASSERT(current_km1 == metrics::km1(_hg));
        _moves.emplace_back(max_gain_hn, from_part, to_part);

        if (current_km1 <= initial_km1) {
          touched_hns_since_improvement = 0;
          last_accepted_index = static_cast<int>(_moves.size()) - 1;
          _gain_manager.resetDelta();
        }
      }
    }

    ASSERT([&]() {
      ASSERT_THAT_PQS_CONTAIN_CORRECT_HYPERNODES();
      return true;
    }());

    // rollback bad moves
    _gain_manager.rollbackDelta();

    int last_index = static_cast<int>(_moves.size()) - 1;
    while (last_index != last_accepted_index) {
      const HypernodeID hn = _moves[last_index].hn;
      const PartitionID from_part = _moves[last_index].to;
      const PartitionID to_part = _moves[last_index].from;

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

      // moving the node back might change the move targets of neighboring nodes
      reinitNodeWithNeighbors(hn, true); // TODO is there a way around this?
    }

    for (const HypernodeID& marked_hn : marked_hns) {
      _hg.unmark(marked_hn);
      resetMinMaxNeighborsFor(marked_hn);
      updateMinMaxNeighborFor(marked_hn);
    }

    ASSERT([&]() {
      ASSERT_THAT_PQS_CONTAIN_CORRECT_HYPERNODES();
      return true;
    }());
    ASSERT(_qg.isAcyclic(), "Rollback produced a cyclic quotient graph!");
    ASSERT(AdjacencyMatrixQuotientGraph<DFSCycleDetector>(_hg, _context).isAcyclic());
    ASSERT(best_metrics.km1 == metrics::km1(_hg));
    best_metrics.cut = metrics::hyperedgeCut(_hg);
    _improved_imbalance += initial_imbalance - best_metrics.imbalance;
    _improved_km1 += initial_km1 - best_metrics.km1;
    return initial_imbalance > best_metrics.imbalance;
  }

  void reinitNodeWithNeighbors(const HypernodeID hn, bool reset = true) {
    if (_qg_changed) {
      return;
    }
    if (reset) {
      _updated_neighbors.resetUsedEntries();
    }

    updateMinMaxNeighborFor(hn);
    _updated_neighbors.set(hn, true);

    for (const HyperedgeID& he : _hg.incidentEdges(hn)) {
      for (const HypernodeID& pin : _hg.pins(he)) {
        if (_updated_neighbors.get(pin)) {
          continue;
        }
        _updated_neighbors.set(pin, true);
        updateMinMaxNeighborFor(pin);
      }
    }
  }

  bool isOverloaded(const PartitionID part) const {
    return _hg.partWeight(part) > _context.partition.max_part_weights[part];
  }

  void initializeMinMaxNeighbors() {
    for (const HypernodeID& hn : _hg.nodes()) {
      updateMinMaxNeighborFor(hn);
    }
  }

  void resetMinMaxNeighborsFor(const HypernodeID hn) {
    _pq_inst[hn].clear();
    if (_hg.active(hn)) {
      _pq.remove(hn, _hg.partID(hn));
    }
    _min_successor_part[hn] = _hg.partID(hn);
    _max_predecessor_part[hn] = _hg.partID(hn);
  }

  void updateMinMaxNeighborFor(const HypernodeID hn) {
    const PartitionID old_min_successor = _min_successor_part[hn];
    const PartitionID old_max_predecessor = _max_predecessor_part[hn];
    const auto& ordering = _qg.topologicalOrdering();
    const auto& inverse_ordering = _qg.inverseTopologicalOrdering();

    _min_successor_part[hn] = _context.partition.k - 1;
    _max_predecessor_part[hn] = 0;

    for (const HyperedgeID& he : _hg.incidentTailEdges(hn)) {
      for (const HypernodeID& head : _hg.heads(he)) {
        _min_successor_part[hn] = std::min<PartitionID>(ordering[_hg.partID(head)], _min_successor_part[hn]);
      }
    }
    for (const HyperedgeID& he : _hg.incidentHeadEdges(hn)) {
      for (const HypernodeID& tail : _hg.tails(he)) {
        ASSERT(ordering[_hg.partID(tail)] <= ordering[_hg.partID(hn)],
               V(ordering[_hg.partID(tail)]) << V(tail) << V(ordering[_hg.partID(hn)]) << V(hn));
        _max_predecessor_part[hn] = std::max<PartitionID>(ordering[_hg.partID(tail)], _max_predecessor_part[hn]);
      }
    }

    ASSERT(_max_predecessor_part[hn] <= ordering[_hg.partID(hn)]);
    ASSERT(_min_successor_part[hn] >= ordering[_hg.partID(hn)]);

    if (old_min_successor != old_max_predecessor) {
      if (!_hg.marked(hn)) {
        _pq.remove(hn, _hg.partID(hn));
      }
      _pq_inst[hn].clear();
    }
    if (_min_successor_part[hn] != _max_predecessor_part[hn]) {
      for (PartitionID topo_part = _max_predecessor_part[hn]; topo_part <= _min_successor_part[hn]; ++topo_part) {
        const PartitionID part = inverse_ordering[topo_part];
        if (part != _hg.partID(hn)) {
          _pq_inst[hn].push(part, _gain_manager.gain(hn, part));
        }
      }
      if (!_hg.marked(hn)) {
        _pq.insert(hn, _hg.partID(hn), _pq_inst[hn].topKey());
      }
    }
  }

  void updateMinMaxNeighborsAfterMove(const HypernodeID moved_hn, const PartitionID from_part,
                                      const PartitionID to_part) {
    _updated_neighbors.resetUsedEntries();
    _min_successor_part[moved_hn] = from_part;
    _max_predecessor_part[moved_hn] = from_part;
    updateMinMaxNeighborFor(moved_hn);
    _updated_neighbors.set(moved_hn, true);


    for (const HyperedgeID& he : _hg.incidentEdges(moved_hn)) {
      for (const HypernodeID& pin : _hg.pins(he)) {
        if (_updated_neighbors.get(pin) || _hg.marked(pin)) {
          continue;
        }
        _updated_neighbors.set(pin, true);
        updateMinMaxNeighborFor(pin);
      }
    }
  }

  void refreshTopologicalOrdering() {
    LOG << "Re-initialize PQs in AcyclicSoftRebalanceRefiner";
    _pq.clear();
    for (const HypernodeID& hn : _hg.nodes()) {
      resetMinMaxNeighborsFor(hn);
    }
    initializeMinMaxNeighbors();
  }

  void move(const HypernodeID moved_hn, const PartitionID from_part, const PartitionID to_part) {
    DBG << "move(" << V(moved_hn) << V(from_part) << V(to_part) << ")";
    ASSERT(_hg.marked(moved_hn));
    ASSERT(_hg.partID(moved_hn) == to_part);

    _updated_neighbors.resetUsedEntries();
    updateMinMaxNeighborsAfterMove(moved_hn, from_part, to_part);

    updateHypernodeGains(moved_hn);
    _updated_neighbors.set(moved_hn, true);

    for (const HyperedgeID& he : _hg.incidentEdges(moved_hn)) {
      for (const HypernodeID& pin : _hg.pins(he)) {
        if (!_updated_neighbors.get(pin)) {
          updateHypernodeGains(pin);
          _updated_neighbors.set(pin, true);
        }
      }
    }

    ASSERT([&]() {
      ASSERT_THAT_PQS_CONTAIN_CORRECT_HYPERNODES();
      return true;
    }());
  }

  void updateHypernodeGains(const HypernodeID hn) {
    if (_hg.marked(hn)) {
      return;
    }
    ASSERT(_hg.active(hn));

    const auto& ordering = _qg.topologicalOrdering();
    const auto& inverse_ordering = _qg.inverseTopologicalOrdering();

    for (PartitionID part_id = ordering[_max_predecessor_part[hn]]; part_id <= ordering[_min_successor_part[hn]]; ++part_id) {
      const PartitionID part = inverse_ordering[part_id];
      if (part != _hg.partID(hn)) {

      }

    }
  }

  bool isMovableTo(const HypernodeID hn, const PartitionID part) const {
    const auto& ordering = _qg.topologicalOrdering();
    return _hg.partID(hn) != part
           && _max_predecessor_part[hn] <= ordering[part]
           && ordering[part] <= _min_successor_part[hn];
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

  void ASSERT_THAT_PQS_CONTAIN_CORRECT_HYPERNODES() const {
    const auto &inverse_ordering = _qg.inverseTopologicalOrdering();

    for (const HypernodeID& hn : _hg.nodes()) {
      if (_max_predecessor_part[hn] != _min_successor_part[hn]) {
        ASSERT(!_pq_inst[hn].empty());
        if (!_hg.marked(hn)) {
          ASSERT(_pq.contains(hn, _hg.partID(hn)));
        }

        std::size_t expected_size = 0;
        HyperedgeWeight max_gain = std::numeric_limits<HyperedgeWeight>::min();

        for (PartitionID topo_part = _max_predecessor_part[hn]; topo_part <= _min_successor_part[hn]; ++topo_part) {
          const PartitionID part = inverse_ordering[topo_part];
          if (part != _hg.partID(hn)) {
            ++expected_size;
            ASSERT(_pq_inst[hn].contains(part));
            max_gain = std::max(max_gain, _pq_inst[hn].getKey(part));
          } else {
            ASSERT(!_pq_inst[hn].contains(part));
          }
        }
        ASSERT(_pq_inst[hn].size() == expected_size);
        ASSERT(max_gain == _pq_inst[hn].topKey());

        if (!_hg.marked(hn)) {
          ASSERT(_pq.contains(hn, _hg.partID(hn)));
          ASSERT(_pq.key(hn, _hg.partID(hn)) == max_gain);
        }
      } else {
        ASSERT(!_pq.contains(hn));
        ASSERT(_pq_inst[hn].empty());
      }
    }
  }
#endif // KAHYPAR_USE_ASSERTIONS

  Hypergraph& _hg;
  const Context& _context;
  AdjacencyMatrixQuotientGraph<DFSCycleDetector>& _qg;
  KMinusOneGainManager& _gain_manager;

  KWayRefinementPQ _pq;
  std::vector<ds::BinaryMaxHeap<PartitionID, Gain>> _pq_inst;

  std::vector<PartitionID> _min_successor_part;
  std::vector<PartitionID> _max_predecessor_part;

  std::vector<Move> _moves;

  ds::FastResetArray<bool> _updated_neighbors;
  bool _qg_changed{false};

  std::size_t _num_refreshes{0};
  std::size_t _num_no_refreshes{0};
  std::size_t _num_moves{0};
  std::size_t _num_positive_gain_moves{0};
  std::size_t _num_negative_gain_moves{0};
  std::size_t _num_zero_gain_moves{0};
  std::size_t _num_moves_in_last_iteration{0};
  HyperedgeWeight _improved_km1{0};
  double _improved_imbalance{0.0};
};
} // namespace kahypar