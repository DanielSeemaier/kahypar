#pragma once

#include <limits>
#include <stack>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "gtest/gtest_prod.h"

#include "kahypar/datastructure/fast_reset_array.h"
#include "kahypar/datastructure/fast_reset_flag_array.h"
#include "kahypar/datastructure/sparse_map.h"
#include "kahypar/datastructure/kway_priority_queue.h"
#include "kahypar/datastructure/binary_heap.h"
#include "kahypar/dag/quotient_graph.h"
#include "kahypar/definitions.h"
#include "kahypar/meta/mandatory.h"
#include "kahypar/meta/template_parameter_to_string.h"
#include "kahypar/partition/context.h"
#include "kahypar/partition/metrics.h"
#include "kahypar/partition/refinement/fm_refiner_base.h"
#include "kahypar/partition/refinement/i_refiner.h"
#include "kahypar/partition/refinement/kway_fm_gain_cache.h"
#include "kahypar/partition/refinement/policies/fm_improvement_policy.h"
#include "kahypar/utils/float_compare.h"
#include "kahypar/utils/randomize.h"
#include "kahypar/utils/timer.h"

namespace kahypar {
template<class StoppingPolicy = Mandatory,
  class FMImprovementPolicy = CutDecreasedOrInfeasibleImbalanceDecreased>
class AcyclicKWayBalanceRefiner final : public IRefiner {
 private:
  static constexpr bool debug = true;

  static constexpr HypernodeID kInvalidHN = std::numeric_limits<HypernodeID>::max();
  static constexpr Gain kInvalidGain = std::numeric_limits<Gain>::min();

  struct PinState {
    char one_pin_in_from_part_before : 1;
    char one_pin_in_to_part_after : 1;
    char two_pins_in_from_part_before : 1;
    char two_pins_in_to_part_after : 1;

    PinState(const bool one_in_from_before, const bool one_in_to_after,
             const bool two_in_from_before, const bool two_in_to_after) :
      one_pin_in_from_part_before(one_in_from_before),
      one_pin_in_to_part_after(one_in_to_after),
      two_pins_in_from_part_before(two_in_from_before),
      two_pins_in_to_part_after(two_in_to_after) {}
  };

#ifdef USE_BUCKET_QUEUE
  using KWayRefinementPQ = ds::KWayPriorityQueue<HypernodeID, Gain,
                                                 std::numeric_limits<Gain>,
                                                 false,
                                                 ds::EnhancedBucketQueue<HypernodeID,
                                                                         Gain,
                                                                         std::numeric_limits<Gain>
                                                                         > >;
#else
  using KWayRefinementPQ = ds::KWayPriorityQueue<HypernodeID, Gain, std::numeric_limits<Gain>>;
#endif

 public:
  AcyclicKWayBalanceRefiner(Hypergraph& hypergraph, const Context& context) :
    _hg(hypergraph),
    _context(context),
    _pq(context.partition.k),
    _tmp_gains(context.partition.k, 0) {}

  ~AcyclicKWayBalanceRefiner() override = default;

  AcyclicKWayBalanceRefiner(const AcyclicKWayBalanceRefiner&) = delete;

  AcyclicKWayBalanceRefiner& operator=(const AcyclicKWayBalanceRefiner&) = delete;

  AcyclicKWayBalanceRefiner(AcyclicKWayBalanceRefiner&&) = delete;

  AcyclicKWayBalanceRefiner& operator=(AcyclicKWayBalanceRefiner&&) = delete;

 private:
  void initializeImpl(const HyperedgeWeight max_gain) override final {
    if (!_is_initialized) {
#ifdef USE_BUCKET_QUEUE
      _pq.initialize(_hg.initialNumNodes(), max_gain);
#else
      unused(max_gain);
      _pq.initialize(_hg.initialNumNodes());
#endif
      _is_initialized = true;
    }
  }

  void performMovesAndUpdateCacheImpl(const std::vector<Move>& moves,
                                      std::vector<HypernodeID>& refinement_nodes,
                                      const UncontractionGainChanges& changes) override final {
  }

  bool refineImpl(std::vector<HypernodeID>& refinement_nodes,
                  const std::array<HypernodeWeight, 2>&,
                  const UncontractionGainChanges&,
                  Metrics& best_metrics) override final {
    _pq.clear();
    _pq_inst.clear();
    _pq_inst.reserve(_hg.initialNumNodes());
    for (HypernodeID hn = 0; hn < _hg.initialNumNodes(); ++hn) {
      _pq_inst.emplace_back(_context.partition.k, 0);
    }
    _hg.resetHypernodeState();

    Randomize::instance().shuffleVector(refinement_nodes, refinement_nodes.size());
    for (const HypernodeID& hn : refinement_nodes) {
      activate(hn);
    }

    // Activate all adjacent free vertices of a fixed vertex in refinement_nodes
    for (const HypernodeID& hn : refinement_nodes) {
      if (_hg.isFixedVertex(hn)) {
        for (const HyperedgeID& he : _hg.incidentEdges(hn)) {
          for (const HypernodeID& pin : _hg.pins(he)) {
            if (!_hg.isFixedVertex(pin) && !_hg.active(pin)) {
              activate(pin);
            }
          }
        }
      }
    }

    const double initial_imbalance = best_metrics.imbalance;
    double current_imbalance;

    QuotientGraph<DFSCycleDetector> qg(_hg, _context);

    BinaryMaxHeap<PartitionID, HypernodeWeight> block_sizes(_context.partition.k, 0);
    for (PartitionID k = 0; k < _context.partition.k; ++k) {
      block_sizes.push(k, _hg.partWeight(k));
    }

    std::size_t num_performed_moves = 0;
    while (!_pq.empty()) {
      Gain max_gain = kInvalidGain;
      HypernodeID max_gain_node = kInvalidHN;
      PartitionID from_part = block_sizes.top();
      while (_pq.empty(from_part) || !_pq.isEnabled(from_part)) {
        LOG << V(from_part) << ": " << V(_pq.empty(from_part)) << V(!_pq.isEnabled(from_part));
        block_sizes.pop();
        if (block_sizes.empty()) {
          break;
        }
        from_part = block_sizes.top();
      }
      if (_pq.empty(from_part) || !_pq.isEnabled(from_part) || _hg.partWeight(from_part) < _context.partition.max_part_weights[from_part]) {
        LOG << _pq.empty(from_part) << " " << !_pq.isEnabled(from_part) << " " << (_hg.partWeight(from_part) < _context.partition.max_part_weights[from_part]);
        break;
      }
      _pq.deleteMaxFromPartition(max_gain_node, max_gain, from_part);
      const PartitionID to_part = _pq_inst[max_gain_node].top();

      ASSERT(_pq_inst[max_gain_node].topKey() == max_gain);
      _pq_inst[max_gain_node].pop();
      ASSERT(!_hg.marked(max_gain_node), V(max_gain_node));
      ASSERT(_hg.isBorderNode(max_gain_node), V(max_gain_node));

      bool do_move = (_hg.partWeight(to_part) + _hg.nodeWeight(max_gain_node)
                      <= _hg.partWeight(from_part) - _hg.nodeWeight(max_gain_node)) && (_hg.partSize(from_part) - 1 != 0);
      bool balance_ok = do_move;
      if (do_move) {
        HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
        do_move = do_move && qg.update(max_gain_node, from_part, to_part);
        HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
        Timer::instance().add(_context, Timepoint::cycle_detector, std::chrono::duration<double>(end - start).count());
      }

      //DBG << V(from_part) << V(to_part) << V(max_gain_node) << V(max_gain) << V(balance_ok) << V(do_move);


      if (do_move) {
        if (block_sizes.contains(from_part)) {
          block_sizes.decreaseKeyBy(from_part, _hg.nodeWeight(max_gain_node));
        }
        block_sizes.increaseKeyBy(to_part, _hg.nodeWeight(max_gain_node));

        _hg.mark(max_gain_node);
        _hg.changeNodePart(max_gain_node, from_part, to_part);
        current_imbalance = metrics::imbalance(_hg, _context);
        ++num_performed_moves;
        LOG << "performed MOVE:" << V(max_gain_node) << V(from_part) << V(to_part) << V(current_imbalance)
            << V(initial_imbalance);

        updateNeighbours(max_gain_node, from_part, to_part);
      } else {
        if (!_pq_inst[max_gain_node].empty()) {
          _pq.insert(max_gain_node, _hg.partID(max_gain_node), _pq_inst[max_gain_node].topKey());
          _pq.enablePart(_hg.partID(max_gain_node));
        }
      }
    }

    DBG << "KWayFM performed" << num_performed_moves << "local search movements";
    ASSERT(qg.isAcyclic(), "Refinement produced a cyclic quotient graph!");
    ASSERT(QuotientGraph<DFSCycleDetector>(_hg, _context).isAcyclic(),
           "Refinement produced a cyclic quotient graph not detected by the quotient graph!");
    return num_performed_moves > 0;
  }

  void updatePin(HypernodeID pin, PartitionID part, HyperedgeWeight delta) {
    ASSERT(_hg.active(pin));
    if (_pq_inst[pin].contains(part)) {
      _pq_inst[pin].increaseKeyBy(part, delta);
    }
    _pq.updateKey(pin, _hg.partID(pin), _pq_inst[pin].topKey());
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void deltaGainUpdates(const HypernodeID pin,
                                                        const PartitionID from_part,
                                                        const PartitionID to_part,
                                                        const HyperedgeID he,
                                                        const HyperedgeWeight he_weight,
                                                        const PinState pin_state) {
    const PartitionID source_part = _hg.partID(pin);
    if (source_part == from_part) {
      if (pin_state.two_pins_in_from_part_before) {
        for (const PartitionID& part : _hg.connectivitySet(he)) {
          updatePin(pin, part, he_weight);
        }
      }
    } else if (source_part == to_part && pin_state.two_pins_in_to_part_after) {
      for (const PartitionID& part : _hg.connectivitySet(he)) {
        updatePin(pin, part, -he_weight);
      }
    }

    if (pin_state.one_pin_in_from_part_before) {
      updatePin(pin, from_part, -he_weight);
    }

    if (pin_state.one_pin_in_to_part_after) {
      updatePin(pin, to_part, he_weight);
    }
  }

  void fullUpdate(const HypernodeID moved_hn, const PartitionID from_part,
                  const PartitionID to_part, const HyperedgeID he) {
    ONLYDEBUG(moved_hn);
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
      if (!_hg.marked(pin)) {
        if (!_hg.active(pin)) {
          _hns_to_activate.push_back(pin);
        } else {
          if (!_hg.isBorderNode(pin)) {
            for (PartitionID k = 0; k < _context.partition.k; ++k) {
              if (_pq.contains(pin, k)) {
                _pq.remove(pin, k);
              }
            }
            _pq_inst[pin].clear();
            _hg.deactivate(pin);
          } else {
            //deltaGainUpdates(pin, from_part, to_part, he, he_weight, pin_state);
            for (PartitionID k = 0; k < _context.partition.k; ++k) {
              if (_pq.contains(pin, k)) {
                _pq.remove(pin, k);
              }
            }
            _pq_inst[pin].clear();
            _hg.deactivate(pin);
            activate(pin);
          }
        }
      }
    }
  }

  void updateNeighbours(const HypernodeID moved_hn, const PartitionID from_part, const PartitionID to_part) {
    for (const HyperedgeID& he : _hg.incidentEdges(moved_hn)) {
      fullUpdate(moved_hn, from_part, to_part, he);
    }
    for (const HypernodeID& hn : _hns_to_activate) {
      if (!_hg.active(hn) && likely(!_hg.isFixedVertex(hn))) {
        activate(hn);
      }
    }
    _hns_to_activate.clear();
  }


  void activate(const HypernodeID hn) {
    if (_hg.isBorderNode(hn) && likely(!_hg.isFixedVertex(hn))) {
      ASSERT(!_hg.active(hn), V(hn));
      ASSERT(!_hg.marked(hn), "Hypernode" << hn << "is already marked");
      insertHNintoPQ(hn);
      _hg.activate(hn);
    }
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void insertHNintoPQ(const HypernodeID hn) {
    ASSERT(_hg.isBorderNode(hn));
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
        continue;
      }
      _pq_inst[hn].push(target_part.key, target_part.value - internal);
    }
    _pq.insert(hn, _hg.partID(hn), _pq_inst[hn].topKey());
    _pq.enablePart(_hg.partID(hn));
  }

  ds::SparseMap<PartitionID, Gain> _tmp_gains;
  std::vector<BinaryMaxHeap<PartitionID, HyperedgeWeight>> _pq_inst;
  KWayRefinementPQ _pq;
  const Context& _context;
  Hypergraph& _hg;
  std::vector<HypernodeID> _hns_to_activate;
};
}  // namespace kahypar
