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
#include "kahypar/partition/refinement/acyclic_km1_refiner.h"

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
  using GainCache = KwayGainCache<Gain>;

  AcyclicKWayBalanceRefiner(Hypergraph& hypergraph, const Context& context) :
    _hg(hypergraph),
    _context(context),
    //_gain_cache(_hg.initialNumNodes(), _context.partition.k),
    _pqs() {}
    //_tmp_gains(_context.partition.k, 0) {}

  ~AcyclicKWayBalanceRefiner() override = default;

  AcyclicKWayBalanceRefiner(const AcyclicKWayBalanceRefiner&) = delete;
  AcyclicKWayBalanceRefiner& operator=(const AcyclicKWayBalanceRefiner&) = delete;

  AcyclicKWayBalanceRefiner(AcyclicKWayBalanceRefiner&&) = delete;
  AcyclicKWayBalanceRefiner& operator=(AcyclicKWayBalanceRefiner&&) = delete;

 private:
  void initializeImpl(const HyperedgeWeight max_gain) override final {
    _is_initialized = true;
//    _gain_cache.clear();
//    for (const HypernodeID& hn : _hg.nodes()) {
//      _tmp_gains.clear();
//      const PartitionID source_part = _hg.partID(hn);
//      HyperedgeWeight internal = 0;
//      for (const HyperedgeID& he : _hg.incidentEdges(hn)) {
//        const HyperedgeWeight he_weight = _hg.edgeWeight(he);
//        internal += _hg.pinCountInPart(he, source_part) != 1 ? he_weight : 0;
//        for (const PartitionID& part : _hg.connectivitySet(he)) {
//          _tmp_gains[part] += he_weight;
//        }
//      }
//
//      for (const auto& target_part : _tmp_gains) {
//        if (target_part.key == source_part) {
//          ASSERT(!_gain_cache.entryExists(hn, source_part), V(hn) << V(source_part));
//          continue;
//        }
//        ASSERT(target_part.value - internal == gainInducedByHypergraph(hn, target_part.key),
//               V(gainInducedByHypergraph(hn, target_part.key)) << V(target_part.value - internal));
//        _gain_cache.initializeEntry(hn, target_part.key, target_part.value - internal);
//      }
//    }
  }

  void performMovesAndUpdateCacheImpl(const std::vector<Move>& moves,
                                      std::vector<HypernodeID>& refinement_nodes,
                                      const UncontractionGainChanges& changes) override final {
  }

  void init(std::vector<HypernodeID>& nodes) {
    _pqs.clear();
    for (PartitionID part = 0; part < _context.partition.k; ++part) {
      _pqs.emplace_back(_context.partition.k);
      _pqs.back().initialize(_hg.initialNumNodes());
    }
    _hg.resetHypernodeState();
  }

  bool refineImpl(std::vector<HypernodeID>& refinement_nodes,
                  const std::array<HypernodeWeight, 2>&,
                  const UncontractionGainChanges&,
                  Metrics& best_metrics) override final {
    double imbalance = best_metrics.imbalance;
    LOG << V(imbalance);

    while (imbalance > _context.partition.epsilon) {
      Randomize::instance().shuffleVector(refinement_nodes, refinement_nodes.size());
      init(refinement_nodes);

      HashMapQuotientGraph<DFSCycleDetector> qg(_hg, _context);
      qg.addMissingEdges();
      qg.reduceToOneRoot();
      const auto qmg = qg.computeMoveGraph();
      const auto ordering = qg.computeStrictTopologicalOrdering();

      for (const HypernodeID& hn : refinement_nodes) {
        std::vector<EdgeWeight> gains(_hg.k());
        const PartitionID source_part = _hg.partID(hn);
        HyperedgeWeight internal = 0;
        for (const HyperedgeID& he : _hg.incidentEdges(hn)) {
          const HyperedgeWeight he_weight = _hg.edgeWeight(he);
          internal += _hg.pinCountInPart(he, source_part) != 1 ? he_weight : 0;
          for (const PartitionID& part : _hg.connectivitySet(he)) {
            gains[part] += he_weight;
          }
        }

        bool activate = false;
        for (const QNodeID& part : qmg.outs(source_part)) {
          if (movable_to(hn, ordering[source_part] < ordering[part])) {
            _pqs[source_part].insert(hn, part, gains[part] - internal);
            _pqs[source_part].enablePart(part);
            activate = true;
          }
        }
        if (activate) {
          _hg.activate(hn);
        }
      }

      std::vector<PartitionID> overloaded_blocks;
      std::vector<PartitionID> underloaded_blocks;
      for (PartitionID k = 0; k < _context.partition.k; ++k) {
        if (_hg.partWeight(k) > _context.partition.max_part_weights[k]) {
          overloaded_blocks.push_back(k);
        } else if (_hg.partWeight(k) <
                   _context.partition.perfect_balance_part_weights[k] * (1.0 - _context.partition.epsilon)) {
          underloaded_blocks.push_back(k);
        }
      }
      if (underloaded_blocks.empty()) {
        for (PartitionID k = 0; k < _context.partition.k; ++k) {
          if (_hg.partWeight(k) < _context.partition.perfect_balance_part_weights[k]) {
            underloaded_blocks.push_back(k);
          }
        }
        LOG << V(overloaded_blocks) << V(underloaded_blocks);
      }
      if (overloaded_blocks.empty() || underloaded_blocks.empty()) {
        LOG << V(overloaded_blocks) << V(underloaded_blocks);
        break;
      }

      std::pair<PartitionID, PartitionID> max_path;
      HyperedgeWeight max_path_gain = std::numeric_limits<HyperedgeWeight>::min();

      for (const PartitionID& from : overloaded_blocks) {
        for (const PartitionID& to : underloaded_blocks) {
          std::vector<QNodeID> path;
          const bool found_path = qmg.findPath(from, to, path);
          ASSERT(found_path, "No path from" << from << "to" << to);
          HyperedgeWeight gain = 0;
          for (std::size_t i = path.size() - 1; i > 0; --i) {
            const QNodeID u = path[i - 1];
            const QNodeID v = path[i];
            if (_pqs[u].empty(v)) { // TODO why.......
              gain = std::numeric_limits<HyperedgeWeight>::min();
              break;
            }
            gain += _pqs[u].maxKey(v);
          }
          if (gain > max_path_gain) {
            max_path_gain = gain;
            max_path = {from, to};
          }
        }
      }

      if (max_path_gain == std::numeric_limits<HyperedgeWeight>::min()) {
        break; // TODO meh..........
      }
      LOG << "Max path from" << max_path.first << "to" << max_path.second << "with gain" << max_path_gain;

      std::vector<QNodeID> path;
      const bool found_path = qmg.findPath(max_path.first, max_path.second, path);
      ASSERT(found_path);
      for (std::size_t i = path.size() - 1; i > 0; --i) {
        const QNodeID u = path[i - 1];
        const QNodeID v = path[i];
        bool move_ok;
        HypernodeID hn;
        HyperedgeWeight gain;
        do { // TODO meh....
          if (_pqs[u].empty(v)) {
            return false; // TODO meh....
          }
          _pqs[u].deleteMaxFromPartition(hn, gain, v);
          move_ok = qg.update(hn, u, v);
          ASSERT(move_ok, "bad move" << V(hn) << V(u) << V(v) << V(_hg.partID(hn))); // TODO not meh....?
        } while (!move_ok);

        ASSERT(move_ok, "bad move" << V(hn) << V(u) << V(v) << V(_hg.partID(hn)));
        _hg.changeNodePart(hn, u, v);
      }

      imbalance = metrics::imbalance(_hg, _context);
      auto km1 = metrics::km1(_hg);
      LOG << V(imbalance) << V(km1) << V(overloaded_blocks) << V(underloaded_blocks);
    }

    return false;
  }

  void fullUpdate(const HypernodeID moved_hn,
                  const PartitionID from_part, const PartitionID to_part,
                  const HyperedgeID he) {
    for (const HypernodeID& pin : _hg.pins(he)) {
      if (!_hg.marked(pin)) {
        if (!_hg.active(pin)) {
          _hns_to_activate.push_back(pin);
        } else {
          for (PartitionID k = 0; k < _context.partition.k; ++k) {
            if (_pqs[from_part].contains(pin, k)) {
              _pqs[from_part].remove(pin, k);
            }
          }
          _hg.deactivate(pin);
          activate(pin);
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

  bool movable_to(const HypernodeID hn, const bool lower) {
    if (lower) {
      for (const HyperedgeID& he : _hg.incidentTailEdges(hn)) {
        for (const HypernodeID& ht : _hg.heads(he)) {
          if (_hg.partID(ht) == _hg.partID(hn)) {
            return false;
          }
        }
      }
    } else {
      for (const HyperedgeID& he : _hg.incidentHeadEdges(hn)) {
        for (const HypernodeID& hh : _hg.tails(he)) {
          if (_hg.partID(hh) == _hg.partID(hn)) {
            return false;
          }
        }
      }
    }
    return true;
  }

  void activate(const HypernodeID hn) {
    std::vector<EdgeWeight> gains(_hg.k());
    const PartitionID source_part = _hg.partID(hn);
    HyperedgeWeight internal = 0;
    for (const HyperedgeID& he : _hg.incidentEdges(hn)) {
      const HyperedgeWeight he_weight = _hg.edgeWeight(he);
      internal += _hg.pinCountInPart(he, source_part) != 1 ? he_weight : 0;
      for (const PartitionID& part : _hg.connectivitySet(he)) {
        gains[part] += he_weight;
      }
    }

    for (PartitionID part = 0; part < _hg.k(); ++part) {
      _pqs[source_part].insert(hn, part, gains[part] - internal);
      _pqs[source_part].enablePart(part);
    }
    _hg.activate(hn);
  }

  std::vector<KWayRefinementPQ> _pqs;
  const Context& _context;
  Hypergraph& _hg;
  std::vector<HypernodeID> _hns_to_activate;
  //GainCache _gain_cache;
  //ds::SparseMap<PartitionID, Gain> _tmp_gains;
};
}  // namespace kahypar
