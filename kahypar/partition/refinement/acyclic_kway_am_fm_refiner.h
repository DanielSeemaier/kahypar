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
template<typename StoppingPolicy = Mandatory,
    typename FMImprovementPolicy = CutDecreasedOrInfeasibleImbalanceDecreased>
class AcyclicKWayAdvancedMovesFMRefiner final : public IRefiner {
private:
  static constexpr bool debug = false;
  static constexpr HypernodeID hn_to_debug = 1709;

public:
  using KWayRefinementPQ = ds::KWayPriorityQueue<HypernodeID, Gain, std::numeric_limits<Gain>, true>;

  AcyclicKWayAdvancedMovesFMRefiner(Hypergraph &hypergraph,
                                    const Context &context,
                                    AdjacencyMatrixQuotientGraph<DFSCycleDetector> &qg,
                                    KMinusOneGainManager &gain_manager) :
      _hg(hypergraph),
      _context(context),
      _qg(qg),
      _gain_manager(gain_manager),
      _stopping_policy(),
      _pq(context.partition.k),
      _min_successor(_hg.initialNumNodes()),
      _max_predecessor(_hg.initialNumNodes()),
      _updated_neighbors(_hg.initialNumNodes(), false) {}

  ~AcyclicKWayAdvancedMovesFMRefiner() override = default;

  AcyclicKWayAdvancedMovesFMRefiner(const AcyclicKWayAdvancedMovesFMRefiner &) = delete;

  AcyclicKWayAdvancedMovesFMRefiner &operator=(const AcyclicKWayAdvancedMovesFMRefiner &) = delete;

  AcyclicKWayAdvancedMovesFMRefiner(AcyclicKWayAdvancedMovesFMRefiner &&) = delete;

  AcyclicKWayAdvancedMovesFMRefiner &operator=(AcyclicKWayAdvancedMovesFMRefiner &&) = delete;

  const std::vector<Move> &moves() const {
    return _moves;
  }

  void printSummary() const override {
    LOG << "[AcyclicKWayAdvancedMovesFMRefiner] Iterations with refresh:" << _num_refreshes;
    LOG << "[AcyclicKWayAdvancedMovesFMRefiner] Iterations without refresh:" << _num_no_refreshes;
    LOG << "[AcyclicKWayAdvancedMovesFMRefiner] Number of moves:" << _num_moves;
    LOG << "[AcyclicKWayAdvancedMovesFMRefiner] Number of moves in the last iteration:" << _num_moves_in_last_iteration;
    LOG << "[AcyclicKWayAdvancedMovesFMRefiner] Positive gain moves:" << _num_positive_gain_moves;
    LOG << "[AcyclicKWayAdvancedMovesFMRefiner] Zero gain moves:" << _num_zero_gain_moves;
    LOG << "[AcyclicKWayAdvancedMovesFMRefiner] Negative gain moves:" << _num_negative_gain_moves;
    LOG << "[AcyclicKWayAdvancedMovesFMRefiner] Improved imbalance by:" << _improved_imbalance;
    LOG << "[AcyclicKWayAdvancedMovesFMRefiner] Improved KM1 by:" << _improved_km1;
  }

  void preUncontraction(const HypernodeID representant) override {
  }

  void postUncontraction(const HypernodeID representant, const std::vector<HypernodeID> &&partners) override {
    if (_qg_changed) { // will recompute anyways
      return;
    }

    ASSERT(partners.size() == 1, "Currently only supports pair contractions.");
    const HypernodeID partner = partners.front();

    for (const HypernodeID hn : {representant, partner}) {
      reinitMaxPredecessorMinSuccessorFor(hn);
      for (const HyperedgeID &he : _hg.incidentTailEdges(hn)) {
        for (const HypernodeID &head : _hg.heads(he)) {
          reinitMaxPredecessorMinSuccessorFor(head);
        }
      }
      for (const HyperedgeID &he : _hg.incidentHeadEdges(hn)) {
        for (const HypernodeID &tail : _hg.tails(he)) {
          reinitMaxPredecessorMinSuccessorFor(tail);
        }
      }
    }

    ASSERT(ASSERT_THAT_MAX_PREDECESSOR_MIN_SUCCESSOR_IS_CORRECT());
  }

private:
  void initializeImpl(const HyperedgeWeight max_gain) final {
    if (!_is_initialized) {
      _pq.initialize(_hg.initialNumNodes());
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

  void performMovesAndUpdateCacheImpl(const std::vector<Move> &moves,
                                      std::vector<HypernodeID> &refinement_nodes,
                                      const UncontractionGainChanges &changes) final {
  }

  bool refineImpl(std::vector<HypernodeID> &refinement_nodes,
                  const std::array<HypernodeWeight, 2> &,
                  const UncontractionGainChanges &,
                  Metrics &best_metrics) final {
    DBG << "refinement_nodes:" << refinement_nodes;

    ASSERT(best_metrics.imbalance == metrics::imbalance(_hg, _context));
    ASSERT(best_metrics.km1 == metrics::km1(_hg));

    _hg.resetHypernodeState();
    _moves.clear();
    _pq.clear();

    if (_qg_changed) {
      refreshTopologicalOrdering();
      ++_num_refreshes;
    } else {
      ++_num_no_refreshes;
    }

    ASSERT(ASSERT_THAT_MAX_PREDECESSOR_MIN_SUCCESSOR_IS_CORRECT());
    Randomize::instance().shuffleVector(refinement_nodes, refinement_nodes.size());
    for (const HypernodeID &ref_hn : refinement_nodes) {
      activate(ref_hn);
    }

    _stopping_policy.resetStatistics();
    std::size_t touched_hns_since_improvement = 0;
    int last_accepted_index = -1;

    const HyperedgeWeight initial_km1 = best_metrics.km1;
    const double initial_imbalance = best_metrics.imbalance;
    HyperedgeWeight current_km1 = best_metrics.km1;
    double current_imbalance = best_metrics.imbalance;

    const double beta = std::log(_hg.currentNumNodes());
    while (!_pq.empty() && !_stopping_policy.searchShouldStop(touched_hns_since_improvement,
                                                              _context, beta, best_metrics.km1, current_km1)) {
      // Select max gain node for movement
      HypernodeID max_gain_hn = Hypergraph::kInvalidHypernodeID;
      PartitionID to_part = Hypergraph::kInvalidPartition;
      Gain max_gain = std::numeric_limits<Gain>::min();
      _pq.deleteMax(max_gain_hn, max_gain, to_part);

      ASSERT(max_gain_hn != Hypergraph::kInvalidHypernodeID);
      ASSERT(_hg.active(max_gain_hn));
      ASSERT(!_hg.marked(max_gain_hn));
      ASSERT(isMovable(max_gain_hn));

      ASSERT(max_gain != std::numeric_limits<Gain>::min());
      ASSERT(max_gain == _gain_manager.gain(max_gain_hn, to_part));
      ASSERT(max_gain == gainInducedByHypergraph(max_gain_hn, to_part));

      ASSERT(to_part != _hg.partID(max_gain_hn));
      ASSERT(to_part == minSuccessor(max_gain_hn) || to_part == maxPredecessor(max_gain_hn));

      const PartitionID from_part = _hg.partID(max_gain_hn);

      // remove all other movements for this hypernode
      if (to_part != minSuccessor(max_gain_hn) && from_part != minSuccessor(max_gain_hn)) {
        ASSERT(_pq.contains(max_gain_hn, minSuccessor(max_gain_hn)));
        _pq.remove(max_gain_hn, minSuccessor(max_gain_hn));
      } else if (to_part != maxPredecessor(max_gain_hn) && from_part != maxPredecessor(max_gain_hn)) {
        ASSERT(_pq.contains(max_gain_hn, maxPredecessor(max_gain_hn)));
        _pq.remove(max_gain_hn, maxPredecessor(max_gain_hn));
      }
      ASSERT(!_pq.contains(max_gain_hn));

      // check if the move is okay
      _hg.mark(max_gain_hn);
      ++touched_hns_since_improvement;
      const bool weight_ok =
          _hg.partWeight(to_part) + _hg.nodeWeight(max_gain_hn) <= _context.partition.max_part_weights[0];
      const bool size_ok = _hg.partSize(from_part) > 1;
      const bool move_ok = weight_ok && size_ok;

      DBG << "Test[" << move_ok << "] HN" << max_gain_hn << "/" << from_part << "-->" << to_part << "/ Gain" << max_gain
          << "/ Weight OK" << weight_ok << "/ Size OK" << size_ok << "/ PQ Size" << _pq.size();

      if (move_ok) {
        ++_num_moves;
        if (_hg.currentNumNodes() == _hg.initialNumNodes()) {
          ++_num_moves_in_last_iteration;
        }

        if (max_gain > 0) {
          ++_num_positive_gain_moves;
        } else if (max_gain < 0) {
          ++_num_negative_gain_moves;
        } else {
          ++_num_zero_gain_moves;
        }

        // update QG (must always work!)
        const bool qg_update_worked = _qg.testAndUpdateBeforeMovement(max_gain_hn, from_part, to_part);
        ASSERT(qg_update_worked);

        // perform the move
        _hg.changeNodePart(max_gain_hn, from_part, to_part);
        _gain_manager.updateAfterMovement(max_gain_hn, from_part, to_part,
                                          [&](const HypernodeID &hn, const PartitionID &part, const Gain &gain) {},
                                          [&](const HypernodeID &hn, const PartitionID &part, const Gain &delta) {
                                            if (_hg.marked(hn) || !_hg.active(hn) ||
                                                part == Hypergraph::kInvalidPartition) {
                                              return;
                                            }
                                            if (part != minSuccessor(hn) && part != maxPredecessor(hn)) {
                                              return;
                                            }

                                            ASSERT(_pq.contains(hn, part));
                                            _pq.updateKeyBy(hn, part, delta);
                                          },
                                          [&](const HypernodeID &hn, const PartitionID &part, const Gain &) {}
        );

        // update movable state (handles connectivity increases / decreases)
        reinitMaxPredecessorMinSuccessorFor(max_gain_hn);
        updateMaxPredecessorMinSuccessorAfterMove(max_gain_hn, from_part, to_part);

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

        const bool improved_km1_within_balance = (current_imbalance <= _context.partition.epsilon) &&
                                                 (current_km1 < best_metrics.km1 || (current_km1 == best_metrics.km1 &&
                                                                                     Randomize::instance().flipCoin()));
        if (improved_km1_within_balance) {
          DBG << "Improved cut to" << current_km1 << "with imbalance" << current_imbalance;

          best_metrics.km1 = current_km1;
          best_metrics.imbalance = current_imbalance;
          _stopping_policy.resetStatistics();
          last_accepted_index = _moves.size();
          touched_hns_since_improvement = 0;
          _gain_manager.resetDelta();
        }

        _moves.emplace_back(max_gain_hn, from_part, to_part);
      } // end move_ok

      // activate neighbors
      for (const HyperedgeID &he : _hg.incidentEdges(max_gain_hn)) {
        for (const HypernodeID &pin : _hg.pins(he)) {
          if (!_hg.marked(pin) && !_hg.active(pin)) {
            activate(pin);
          }
        }
      }

      ASSERT(ASSERT_THAT_MAX_PREDECESSOR_MIN_SUCCESSOR_IS_CORRECT());
      ASSERT(ASSERT_THAT_CONTAINED_HYPERNODES_ARE_CORRECT());
    }

    // rollback bad moves
    _gain_manager.rollbackDelta();

    int last_index = static_cast<int>(_moves.size()) - 1;
    while (last_index != last_accepted_index) {
      const HypernodeID hn = _moves[last_index].hn;
      const PartitionID from_part = _moves[last_index].to;
      const PartitionID to_part = _moves[last_index].from;
      DBG << "Rollback" << hn << "/" << from_part << "-->" << to_part;

      const bool success = _qg.testAndUpdateBeforeMovement(hn, from_part, to_part);
      ASSERT(success);

      _hg.changeNodePart(hn, from_part, to_part);

      _moves.pop_back();
      --last_index;
      --_num_moves;

      if (_hg.currentNumNodes() == _hg.initialNumNodes()) {
        --_num_moves_in_last_iteration;
      }

      reinitMaxPredecessorMinSuccessorFor(hn);
      updateMaxPredecessorMinSuccessorAfterMove<false>(hn, from_part, to_part);
    }

    ASSERT(ASSERT_THAT_MAX_PREDECESSOR_MIN_SUCCESSOR_IS_CORRECT());
    ASSERT(_qg.isAcyclic(), "Rollback produced a cyclic quotient graph!");
    ASSERT(AdjacencyMatrixQuotientGraph<DFSCycleDetector>(_hg, _context).isAcyclic());
    ASSERT(best_metrics.km1 == metrics::km1(_hg));
    best_metrics.cut = metrics::hyperedgeCut(_hg);
    _improved_imbalance += initial_imbalance - best_metrics.imbalance;
    _improved_km1 += initial_km1 - best_metrics.km1;

    return FMImprovementPolicy::improvementFound(best_metrics.km1, initial_km1,
                                                 best_metrics.imbalance, initial_imbalance,
                                                 _context.partition.epsilon);
  }

  void activate(const HypernodeID hn) {
    if (isMovable(hn)) {
      DBG << "Activating" << hn;
      ASSERT(!_hg.active(hn));
      ASSERT(!_hg.marked(hn));

      if (_hg.partID(hn) != minSuccessor(hn)) {
        insertIntoPQ(hn, minSuccessor(hn));
      }
      if (_hg.partID(hn) != maxPredecessor(hn)) {
        insertIntoPQ(hn, maxPredecessor(hn));
      }
      _hg.activate(hn);
    }
  }

  void insertIntoPQ(const HypernodeID hn, const PartitionID part) {
    _pq.insert(hn, part, _gain_manager.gain(hn, part));
    if (_hg.partWeight(part) < _context.partition.max_part_weights[0]) {
      _pq.enablePart(part);
    }
  }

  bool isFull(const PartitionID part) const {
    return _hg.partWeight(part) >= _context.partition.max_part_weights[0];
  }

  void refreshTopologicalOrdering() {
    for (const HypernodeID &hn : _hg.nodes()) {
      reinitMaxPredecessorMinSuccessorFor(hn);
    }
    ASSERT(ASSERT_THAT_MAX_PREDECESSOR_MIN_SUCCESSOR_IS_CORRECT());
  }

  bool reinitMaxPredecessorFor(const HypernodeID hn) {
    const auto &inv_top = _qg.inverseTopologicalOrdering();
    PartitionID max_predecessor = -1;

    for (const HyperedgeID &he : _hg.incidentHeadEdges(hn)) {
      for (const HypernodeID &tail : _hg.tails(he)) {
        const PartitionID part = inv_top[_hg.partID(tail)];
        max_predecessor = std::max(max_predecessor, part);
      }
    }

    const auto old_max_predecessor = _max_predecessor[hn];
    if (max_predecessor == -1) {
      _max_predecessor[hn] = inv_top[_hg.partID(hn)];
    } else {
      _max_predecessor[hn] = max_predecessor;
    }
    DBGC(hn == hn_to_debug)
    << "reinitMaxPredecessorFor changed" << V(old_max_predecessor) << "to" << V(_max_predecessor[hn]);
    return _max_predecessor[hn] != old_max_predecessor;
  }

  bool reinitMinSuccessorFor(const HypernodeID hn) {
    const auto &inv_top = _qg.inverseTopologicalOrdering();
    PartitionID min_successor = std::numeric_limits<PartitionID>::max();

    for (const HyperedgeID &he : _hg.incidentTailEdges(hn)) {
      for (const HypernodeID &head : _hg.heads(he)) {
        const PartitionID part = inv_top[_hg.partID(head)];
        min_successor = std::min(min_successor, part);
      }
    }

    const auto old_min_successor = _min_successor[hn];
    if (min_successor == std::numeric_limits<PartitionID>::max()) {
      _min_successor[hn] = inv_top[_hg.partID(hn)];
    } else {
      _min_successor[hn] = min_successor;
    }
    DBGC(hn == hn_to_debug) << "reinitMinSuccessorFor changed" << V(old_min_successor) << "to" << V(_min_successor[hn]);
    return old_min_successor != _min_successor[hn];
  }

  bool reinitMaxPredecessorMinSuccessorFor(const HypernodeID &hn) {
    const auto max_predecessor_changed = reinitMaxPredecessorFor(hn);
    const auto min_successor_changed = reinitMinSuccessorFor(hn);
    return max_predecessor_changed || min_successor_changed;
  }

  template<bool update_pq = true>
  void updateMaxPredecessorMinSuccessorAfterMove(const HypernodeID hn, const PartitionID from_part,
                                                 const PartitionID to_part) {
    const auto &inv_top = _qg.inverseTopologicalOrdering();
    const auto &top = _qg.topologicalOrdering();
    const bool forward_move = inv_top[from_part] < inv_top[to_part];

    if (forward_move) {
      // successor
      for (const HyperedgeID &he : _hg.incidentTailEdges(hn)) {
        for (const HypernodeID &head : _hg.heads(he)) {
          DBGC(head == hn_to_debug) << "forward_move, succ," << V(head) << V(_hg.partID(head))
                                    << V(_max_predecessor[head]) << V(maxPredecessor(head));
          ASSERT(inv_top[_hg.partID(head)] >= inv_top[to_part]);
          if (maxPredecessor(head) == to_part) continue; // nothing changes

          const auto old_max_predecessor = maxPredecessor(head);

          if (_hg.partID(head) == to_part) { // moving it to head's block, max pred was != to_part
            _max_predecessor[head] = inv_top[to_part];

            if (update_pq && _hg.active(head)) { // no longer backward movable
              ASSERT(_pq.contains(head, old_max_predecessor));
              _pq.remove(head, old_max_predecessor);
              if (!isMovable(head)) _hg.deactivate(head);
            }
          } else {
            const auto backward_movable_before = _hg.partID(head) != maxPredecessor(head);
            const auto max_predecessor_changed = reinitMaxPredecessorFor(head);
            const auto backward_movable_after = _hg.partID(head) != maxPredecessor(head);
            if (!update_pq || !_hg.active(head)) continue; // deactivated HNs get activated separately

            ASSERT(!(!backward_movable_before && backward_movable_after));

            if (backward_movable_before && !backward_movable_after) {
              ASSERT(_pq.contains(head, old_max_predecessor));
              _pq.remove(head, old_max_predecessor);
              if (!isMovable(head)) _hg.deactivate(head);
            } else if (max_predecessor_changed) {
              ASSERT(backward_movable_before && backward_movable_after);
              ASSERT(_pq.contains(head, old_max_predecessor));
              _pq.remove(head, old_max_predecessor);
              insertIntoPQ(head, maxPredecessor(head));
            }
          }
        }
      }

      // predecessor
      for (const HyperedgeID &he : _hg.incidentHeadEdges(hn)) {
        for (const HypernodeID &tail : _hg.tails(he)) {
          DBGC(tail == hn_to_debug) << "forward_move, pred," << V(tail) << V(_hg.partID(tail))
                                    << V(_min_successor[tail]) << V(minSuccessor(tail));
          ASSERT(inv_top[_hg.partID(tail)] <= inv_top[from_part]);
          if (minSuccessor(tail) != from_part) continue; // nothing changes

          const auto forward_movable_before = _hg.partID(tail) != minSuccessor(tail);
          const auto min_successor_changed = reinitMinSuccessorFor(tail);
          const auto forward_movable_after = _hg.partID(tail) != minSuccessor(tail);
          if (!update_pq || !_hg.active(tail)) continue; // deactivated HNs get activated separately

          ASSERT(!forward_movable_before || forward_movable_after);

          if (!forward_movable_before && forward_movable_after) {
            ASSERT(!_pq.contains(tail, from_part));
            insertIntoPQ(tail, minSuccessor(tail));
          } else if (min_successor_changed) {
            ASSERT(forward_movable_before && forward_movable_after);
            ASSERT(_pq.contains(tail, from_part));
            _pq.remove(tail, from_part);
            insertIntoPQ(tail, minSuccessor(tail));
          }
        }
      }
    } else {
      // successor
      for (const HyperedgeID &he : _hg.incidentTailEdges(hn)) {
        for (const HypernodeID &head : _hg.heads(he)) {
          DBGC(head == hn_to_debug) << "backward_move, succ," << V(head) << V(_hg.partID(head))
                                    << V(_max_predecessor[head]) << V(maxPredecessor(head));
          ASSERT(inv_top[_hg.partID(head)] >= inv_top[from_part]);
          if (_max_predecessor[head] != inv_top[from_part]) continue;

          const auto backward_movable_before = _hg.partID(head) != maxPredecessor(head);
          const auto max_predecessor_changed = reinitMaxPredecessorFor(head);
          const auto backward_movable_after = _hg.partID(head) != maxPredecessor(head);
          if (!update_pq || !_hg.active(head)) continue; // deactivated HNs get activated separately
          DBGC(head == hn_to_debug) << "\t" << V(backward_movable_before) << V(backward_movable_after)
                                    << V(_max_predecessor[head]) << V(maxPredecessor(head));

          ASSERT(!backward_movable_before || backward_movable_after);

          if (!backward_movable_before && backward_movable_after) {
            ASSERT(!_pq.contains(head, maxPredecessor(head)));
            insertIntoPQ(head, maxPredecessor(head));
          } else if (max_predecessor_changed) {
            ASSERT(backward_movable_before && backward_movable_after);
            ASSERT(_pq.contains(head, from_part));
            _pq.remove(head, from_part);
            insertIntoPQ(head, maxPredecessor(head));
          }
        }
      }

      // predecessor
      for (const HyperedgeID &he : _hg.incidentHeadEdges(hn)) {
        for (const HypernodeID &tail : _hg.tails(he)) {
          DBGC(tail == hn_to_debug) << "backward_move, pred," << V(tail) << V(_hg.partID(tail))
                                    << V(_min_successor[tail]) << V(minSuccessor(tail)) << V(to_part);
          ASSERT(inv_top[_hg.partID(tail)] <= inv_top[to_part]);
          if (minSuccessor(tail) == to_part) continue; // nothing changes
          const auto old_min_successor = minSuccessor(tail);

          if (_hg.partID(tail) == to_part) {
            _min_successor[tail] = inv_top[to_part];

            if (update_pq && _hg.active(tail)) {
              ASSERT(_pq.contains(tail, old_min_successor));
              _pq.remove(tail, old_min_successor);
              if (!isMovable(tail)) _hg.deactivate(tail);
            }
          } else {
            const auto forward_movable_before = _hg.partID(tail) != minSuccessor(tail);
            const auto min_successor_changed = reinitMinSuccessorFor(tail);
            const auto forward_movable_after = _hg.partID(tail) != minSuccessor(tail);
            DBGC(tail == hn_to_debug) << "State:" << V(forward_movable_before) << V(min_successor_changed) << V(forward_movable_after);
            if (!update_pq || !_hg.active(tail)) continue;

            if (forward_movable_before && !forward_movable_after) {
              ASSERT(_pq.contains(tail, old_min_successor));
              _pq.remove(tail, old_min_successor);
              if (!isMovable(tail)) _hg.deactivate(tail);
            } else if (min_successor_changed) {
              ASSERT(forward_movable_before && forward_movable_after);
              ASSERT(_pq.contains(tail, old_min_successor));
              _pq.remove(tail, old_min_successor);
              insertIntoPQ(tail, minSuccessor(tail));
            }
          }

        }
      }
    }
  }

  bool isMovable(const HypernodeID hn) const {
    return maxPredecessor(hn) != minSuccessor(hn);
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
    for (const HyperedgeID &he : _hg.incidentEdges(hn)) {
      gain += gainInducedByHyperedge(hn, he, target_part);
    }
    return gain;
  }

  bool isConnectedTo(const HypernodeID pin, const PartitionID part) const {
    for (const HyperedgeID &he : _hg.incidentEdges(pin)) {
      if (_hg.pinCountInPart(he, part) > 0) {
        return true;
      }
    }
    return false;
  }

  bool ASSERT_THAT_CONTAINED_HYPERNODES_ARE_CORRECT() const {
    const auto &topo = _qg.inverseTopologicalOrdering();

    for (const HypernodeID &hn : _hg.nodes()) {
      const auto part = _hg.partID(hn);
      const auto min_successor = minSuccessor(hn);
      const auto max_predecessor = maxPredecessor(hn);

      if (min_successor != max_predecessor) {
        if (_hg.marked(hn) || !_hg.active(hn)) {
          ASSERT(!_pq.contains(hn));
          continue;
        }

        if (min_successor != part) {
          ASSERT(_pq.contains(hn, min_successor),
                 V(hn) << V(min_successor) << V(part) << "Top" << _qg.topologicalOrdering());
          ASSERT(_pq.key(hn, min_successor) == _gain_manager.gain(hn, min_successor));
        }
        if (max_predecessor != part) {
          ASSERT(_pq.contains(hn, max_predecessor));
          ASSERT(_pq.key(hn, max_predecessor) == _gain_manager.gain(hn, max_predecessor));
        }
      }
    }

    return true;
  }

  bool ASSERT_THAT_MAX_PREDECESSOR_MIN_SUCCESSOR_IS_CORRECT() const {
    const auto &inv_top = _qg.inverseTopologicalOrdering();
    const auto &top = _qg.topologicalOrdering();

    for (const HypernodeID &hn : _hg.nodes()) {
      PartitionID max_predecessor = -1;
      PartitionID min_successor = std::numeric_limits<PartitionID>::max();
      for (const HyperedgeID &he : _hg.incidentTailEdges(hn)) {
        for (const HypernodeID &head : _hg.heads(he)) {
          const PartitionID part = inv_top[_hg.partID(head)];
          min_successor = std::min(min_successor, part);
        }
      }
      for (const HyperedgeID &he : _hg.incidentHeadEdges(hn)) {
        for (const HypernodeID &tail : _hg.tails(he)) {
          const PartitionID part = inv_top[_hg.partID(tail)];
          max_predecessor = std::max(max_predecessor, part);
        }
      }

      if (max_predecessor == -1) {
        ASSERT(maxPredecessor(hn) == _hg.partID(hn), hn);
      } else {
        ASSERT(maxPredecessor(hn) == top[max_predecessor],
               V(hn) << V(maxPredecessor(hn)) << V(top[max_predecessor]) << V(_hg.partID(hn)));
      }
      if (min_successor == std::numeric_limits<PartitionID>::max()) {
        ASSERT(minSuccessor(hn) == _hg.partID(hn), hn);
      } else {
        ASSERT(minSuccessor(hn) == top[min_successor],
               V(hn) << V(minSuccessor(hn)) << V(top[min_successor]) << V(_hg.partID(hn)));
      }
    }

    return true;
  }

#endif // KAHYPAR_USE_ASSERTIONS

  PartitionID maxPredecessor(const HypernodeID &hn) const {
    return _qg.topologicalOrdering()[_max_predecessor[hn]];
  }

  PartitionID minSuccessor(const HypernodeID &hn) const {
    return _qg.topologicalOrdering()[_min_successor[hn]];
  }

  Hypergraph &_hg;
  const Context &_context;
  AdjacencyMatrixQuotientGraph<DFSCycleDetector> &_qg;
  KMinusOneGainManager &_gain_manager;
  StoppingPolicy _stopping_policy;

  KWayRefinementPQ _pq;

  std::vector<PartitionID> _max_predecessor;
  std::vector<PartitionID> _min_successor;

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