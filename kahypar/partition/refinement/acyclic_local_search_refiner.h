#pragma once

#include <array>

#include "kahypar/datastructure/fast_reset_array.h"
#include "kahypar/datastructure/fast_reset_flag_array.h"
#include "kahypar/datastructure/sparse_map.h"
#include "kahypar/dag/quotient_graph.h"
#include "kahypar/definitions.h"
#include "kahypar/meta/mandatory.h"
#include "kahypar/partition/refinement/policies/fm_improvement_policy.h"
#include "kahypar/meta/template_parameter_to_string.h"
#include "kahypar/partition/context.h"
#include "kahypar/partition/metrics.h"
#include "kahypar/partition/refinement/i_refiner.h"
#include "kahypar/partition/refinement/km1_gain_manager.h"
#include "kahypar/datastructure/kway_priority_queue.h"
#include "kahypar/utils/randomize.h"
#include "kahypar/utils/timer.h"

namespace kahypar {
template<typename StoppingPolicy = Mandatory,
         typename FMImprovementPolicy = CutDecreasedOrInfeasibleImbalanceDecreased>
class AcyclicLocalSearchRefiner final : public IRefiner {
 private:
  using KWayRefinementPQ = ds::KWayPriorityQueue<HypernodeID, Gain, std::numeric_limits<Gain>>;

  static constexpr bool debug = true;
  static constexpr HypernodeID hn_to_debug = 0;

  static constexpr HypernodeID kInvalidHN = std::numeric_limits<HypernodeID>::max();
  static constexpr Gain kInvalidGain = std::numeric_limits<Gain>::min();
  static constexpr HyperedgeWeight kInvalidDecrease = std::numeric_limits<PartitionID>::min();

 public:
  AcyclicLocalSearchRefiner(Hypergraph& hypergraph, const Context& context,
                              AdjacencyMatrixQuotientGraph<DFSCycleDetector>& qg,
                              KMinusOneGainManager& gain_manager) :
    _hg(hypergraph),
    _context(context),
    _pq(context.partition.k),
    _qg(qg),
    _updated_neighbors(_hg.initialNumNodes(), false),
    _gain_manager(gain_manager),
    _stopping_policy() {}

  ~AcyclicLocalSearchRefiner() override = default;

  AcyclicLocalSearchRefiner(const AcyclicLocalSearchRefiner&) = delete;
  AcyclicLocalSearchRefiner& operator=(const AcyclicLocalSearchRefiner&) = delete;

  AcyclicLocalSearchRefiner(AcyclicLocalSearchRefiner&&) = delete;
  AcyclicLocalSearchRefiner& operator=(AcyclicLocalSearchRefiner&&) = delete;

  const std::vector<Move>& moves() {
    return _moves;
  }

 private:
  void initializeImpl(const HyperedgeWeight max_gain) final {
    if (!_is_initialized) {
      _pq.initialize(_hg.initialNumNodes());
      _is_initialized = true;
    }
  }

  void performMovesAndUpdateCacheImpl(const std::vector<Move>& moves,
                                      std::vector<HypernodeID>& refinement_nodes,
                                      const UncontractionGainChanges& changes) final {
  }

  bool refineImpl(std::vector<HypernodeID>& refinement_nodes,
                  const std::array<HypernodeWeight, 2>&,
                  const UncontractionGainChanges&,
                  Metrics& best_metrics) final {
    ASSERT(best_metrics.km1 == metrics::km1(_hg));
    ASSERT(best_metrics.imbalance == metrics::imbalance(_hg, _context));

    _hg.resetHypernodeState();
    _moves.clear();
    _pq.clear();

    for (const HypernodeID& hn : refinement_nodes) {
      activate(hn);
    }

    ASSERT([&]() {
      ASSERT_THAT_HYPERNODES_CONTAINED_IN_PQ_ARE_CORRECT();
      return true;
    }());

    const double initial_imbalance = best_metrics.imbalance;
    double current_imbalance = best_metrics.imbalance;
    const HyperedgeWeight initial_km1 = best_metrics.km1;
    HyperedgeWeight current_km1 = best_metrics.km1;

    int min_cut_index = -1;
    int touched_hns_since_last_improvement = 0;
    _stopping_policy.resetStatistics();

    const double beta = std::log(_hg.currentNumNodes());
    while (!_pq.empty() && !_stopping_policy.searchShouldStop(touched_hns_since_last_improvement,
                                                              _context, beta, best_metrics.km1, current_km1)) {
      // Step 1: retrieve next movement
      Gain max_gain = kInvalidGain;
      HypernodeID max_gain_node = kInvalidHN;
      PartitionID to_part = Hypergraph::kInvalidPartition;
      _pq.deleteMax(max_gain_node, max_gain, to_part);

      ASSERT(max_gain != kInvalidGain);
      ASSERT(max_gain_node != kInvalidHN);
      ASSERT(to_part != Hypergraph::kInvalidPartition);
      ASSERT(_hg.active(max_gain_node));
      ASSERT(!_hg.marked(max_gain_node));
      ASSERT(max_gain == _gain_manager.adjacentGain(max_gain_node, to_part));
      ASSERT(max_gain == gainInducedByHypergraph(max_gain_node, to_part));
      ASSERT(_hg.isBorderNode(max_gain_node));
      ASSERT(isConnectedTo(max_gain_node, to_part));

      const PartitionID from_part = _hg.partID(max_gain_node);

      // Step 2: remove all movements for this hypernode
      for (const PartitionID& part : _gain_manager.adjacentParts(max_gain_node)) {
        if (part == to_part) {
          ASSERT(!_pq.contains(max_gain_node, part));
          continue;
        }
        ASSERT(_pq.contains(max_gain_node, part));
        _pq.remove(max_gain_node, part);
      }

      // Step 3: check if the move is okay
      _hg.mark(max_gain_node);
      ++touched_hns_since_last_improvement;

      bool move_ok = ((_hg.partWeight(to_part) + _hg.nodeWeight(max_gain_node))
                      <= _context.partition.max_part_weights[to_part]) // to_part does not become overloaded
                     && (_hg.partSize(from_part) > 1); // from_part does not become empty
      const bool weights_ok = move_ok;

      // quotient graph remains acyclic
      if (move_ok) {
        HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
        move_ok = move_ok && _qg.testAndUpdateBeforeMovement(max_gain_node, from_part, to_part);
        HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
        Timer::instance().add(_context, Timepoint::cycle_detector,
                              std::chrono::duration<double>(end - start).count());
      }

      DBG << "Move[" << move_ok << "]" << max_gain_node << "/" << from_part << "-->" << to_part << "/"
          << max_gain << "/ Weights OK:" << weights_ok << "/ Move OK:" << move_ok;

      // Step 4: perform movement
      if (move_ok) {
        _hg.changeNodePart(max_gain_node, from_part, to_part);

        // enable / disable underloaded / overloaded blocks
        if (_hg.partWeight(to_part) >= _context.partition.max_part_weights[to_part]) {
          _pq.disablePart(to_part);
        }
        if (_hg.partWeight(from_part) < _context.partition.max_part_weights[from_part]) {
          _pq.enablePart(from_part);
        }

        current_imbalance = metrics::imbalance(_hg, _context);
        current_km1 -= max_gain;
        ASSERT(current_km1 == metrics::km1(_hg));
        _stopping_policy.updateStatistics(max_gain);

        // perform gain updates
        _hns_to_activate.clear();
        _gain_manager.updateAfterMovement(max_gain_node, from_part, to_part,
                                          [&](const HypernodeID& hn, const PartitionID& part, const Gain& gain) {
                                            if (_hg.marked(hn)) {
                                              return;
                                            }

                                            // either node activation
                                            if (!_hg.active(hn)) {
                                              _hns_to_activate.push_back(hn);
                                              return;
                                            }

                                            // or just connectivity increase
                                            _pq.insert(hn, to_part, gain);
                                            if (_hg.partWeight(to_part) < _context.partition.max_part_weights[0]) {
                                              _pq.enablePart(to_part);
                                            }
                                          },
                                          [&](const HypernodeID& hn, const PartitionID& part, const Gain& delta) {
                                            if (_hg.marked(hn) || !_hg.active(hn) ||
                                                part == Hypergraph::kInvalidPartition) {
                                              return;
                                            }

                                            ASSERT(_pq.contains(hn, part));
                                            _pq.updateKeyBy(hn, part, delta);
                                          },
                                          [&](const HypernodeID& hn, const PartitionID& part, const Gain&) {
                                            if (_hg.marked(hn) || !_hg.active(hn)) {
                                              return;
                                            }

                                            ASSERT(_pq.contains(hn, part));
                                            _pq.remove(hn, part);
                                          });

        for (const HypernodeID& hn_to_activate : _hns_to_activate) {
          if (!_hg.active(hn_to_activate)) {
            activate(hn_to_activate);
          }
        }

        // check if we accept the new state
        const bool improved_km1_within_balance = (current_imbalance <= _context.partition.epsilon) &&
                                                 (current_km1 < best_metrics.km1);
        bool improved_balance_less_equal_km1 = current_imbalance < best_metrics.imbalance;
        if (!_context.imbalanced_intermediate_step) {
          improved_balance_less_equal_km1 &= current_km1 <= best_metrics.km1;
        }

        if (improved_km1_within_balance || improved_balance_less_equal_km1) {
          DBGC(max_gain == 0) << "KWayFM improved balance between" << from_part
                              << "and" << to_part << "(max_gain=" << max_gain << ")";
          DBGC(current_km1 < best_metrics.km1) << "KWayFM improved cut from " << best_metrics.km1
                                               << "to" << current_km1;
          best_metrics.km1 = current_km1;
          best_metrics.imbalance = current_imbalance;
          _stopping_policy.resetStatistics();
          min_cut_index = _moves.size();
          touched_hns_since_last_improvement = 0;
          _gain_manager.resetDelta();
        }
        _moves.emplace_back(max_gain_node, from_part, to_part);

        ASSERT([&]() {
          ASSERT_THAT_HYPERNODES_CONTAINED_IN_PQ_ARE_CORRECT();
          return true;
        }());
      }
    }

    // finally, rollback to the best last accepted state
    int last_index = _moves.size() - 1;
    while (last_index != min_cut_index) {
      const HypernodeID hn = _moves[last_index].hn;
      const PartitionID from_part = _moves[last_index].to;
      const PartitionID to_part = _moves[last_index].from;
      DBG << "Rollback" << hn << "/" << from_part << "-->" << to_part;

      const bool success = _qg.testAndUpdateBeforeMovement(hn, from_part, to_part);
      ASSERT(success);
      _hg.changeNodePart(hn, from_part, to_part);

      _moves.pop_back();
      --last_index;
    }
    _gain_manager.rollbackDelta();

    ASSERT(_qg.isAcyclic());
    ASSERT(AdjacencyMatrixQuotientGraph<DFSCycleDetector>(_hg, _context).isAcyclic());
    ASSERT(best_metrics.km1 == metrics::km1(_hg));
    best_metrics.imbalance = metrics::imbalance(_hg, _context);

    DBG << "Improved KM1:" << initial_km1 << "-->" << best_metrics.km1;
    DBG << "Changed imbalance:" << initial_imbalance << "-->" << best_metrics.imbalance;

    return FMImprovementPolicy::improvementFound(best_metrics.km1, initial_km1,
                                                 best_metrics.imbalance, initial_imbalance,
                                                 _context.partition.epsilon);
  }

  void activate(const HypernodeID hn) {
    if (_hg.isBorderNode(hn)) {
      ASSERT(!_hg.active(hn));
      ASSERT(!_hg.marked(hn));

      insertHypernodeIntoPQ(hn);
      _hg.activate(hn);
    }
  }

  void insertHypernodeIntoPQ(const HypernodeID hn) {
    ASSERT(_hg.isBorderNode(hn));

    for (const PartitionID& part : _gain_manager.adjacentParts(hn)) {
      ASSERT(part != _hg.partID(hn));
      ASSERT(_gain_manager.adjacentGain(hn, part) == gainInducedByHypergraph(hn, part));
      ASSERT(isConnectedTo(hn, part));

      DBG << "Insert HN" << hn << "into" << part;

      _pq.insert(hn, part, _gain_manager.adjacentGain(hn, part));
      if (_hg.partWeight(part) < _context.partition.max_part_weights[0]) {
        _pq.enablePart(part);
      }
    }
  }

#ifdef KAHYPAR_USE_ASSERTIONS
  void ASSERT_THAT_HYPERNODES_CONTAINED_IN_PQ_ARE_CORRECT() const {
    for (const HypernodeID& hn : _hg.nodes()) {
      if (!_hg.isBorderNode(hn) || _hg.marked(hn) || !_hg.active(hn)) {
        ASSERT(!_pq.contains(hn));
        continue;
      }

      ASSERT(_hg.active(hn));
      ASSERT(_pq.contains(hn));

      for (PartitionID part = 0; part < _context.partition.k; ++part) {
        if (part == _hg.partID(hn)) {
          ASSERT(!_pq.contains(hn, part));
          continue;
        }
        if (isConnectedTo(hn, part)) {
          ASSERT(_pq.contains(hn, part), "Active HN not contained in PQ:" << V(hn) << V(_hg.partID(hn)) << V(part));
          ASSERT(_pq.key(hn, part) == _gain_manager.adjacentGain(hn, part));
          ASSERT(_pq.key(hn, part) == gainInducedByHypergraph(hn, part));
        } else {
          ASSERT(!_pq.contains(hn, part));
        }
      }
    }
  }

  Gain gainInducedByHypergraph(const HypernodeID hn, const PartitionID target_part) const {
    Gain gain = 0;
    for (const HyperedgeID& he : _hg.incidentEdges(hn)) {
      // this might fail in test cases where we don't remove singleton hyperegdes
      // ASSERT(_hg.edgeSize(he) > 1, V(he));
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

  bool isConnectedTo(const HypernodeID hn, const PartitionID part) const {
    for (const HyperedgeID& he : _hg.incidentEdges(hn)) {
      if (_hg.pinCountInPart(he, part) > 0) {
        return true;
      }
    }

    return false;
  }
#endif

  Hypergraph& _hg;
  const Context& _context;
  KWayRefinementPQ _pq;
  AdjacencyMatrixQuotientGraph<DFSCycleDetector>& _qg;
  std::vector<Move> _moves;
  ds::FastResetArray<bool> _updated_neighbors;
  KMinusOneGainManager& _gain_manager;
  StoppingPolicy _stopping_policy;
  std::vector<HypernodeID> _hns_to_activate;
};
} // namespace kahypar