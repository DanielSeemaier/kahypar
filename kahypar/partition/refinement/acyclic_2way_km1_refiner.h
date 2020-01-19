#pragma once

#include <array>

#include "kahypar/partition/refinement/policies/fm_improvement_policy.h"
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
#include "kahypar/utils/htimer.h"

namespace kahypar {
class AcyclicTwoWayKMinusOneRefiner final : public IRefiner {
 private:
  using KWayRefinementPQ = ds::KWayPriorityQueue<HypernodeID, Gain, std::numeric_limits<Gain>>;

  static constexpr bool debug = false;
  static constexpr HypernodeID hn_to_debug = 10111;

  static constexpr std::size_t SUCCESSORS = 0;
  static constexpr std::size_t PREDECESSORS = 1;

 public:
  AcyclicTwoWayKMinusOneRefiner(Hypergraph& hypergraph, const Context& context,
                                std::shared_ptr<AdjacencyMatrixQuotientGraph<DFSCycleDetector>> qg,
                                std::shared_ptr<KMinusOneGainManager> gain_manager) :
    _hg(hypergraph),
    _context(context),
    _qg(std::move(qg)),
    _pqs{ds::BinaryMaxHeap<HypernodeID, Gain>(_hg.initialNumNodes()),
         ds::BinaryMaxHeap<HypernodeID, Gain>(_hg.initialNumNodes())},
    _fixtures{std::vector<HypernodeID>(_hg.initialNumNodes()),
              std::vector<HypernodeID>(_hg.initialNumNodes())},
    _gain_manager(std::move(gain_manager)),
    _control_qg_and_gm(false) {}

  AcyclicTwoWayKMinusOneRefiner(Hypergraph& hypergraph, const Context& context) :
    _hg(hypergraph),
    _context(context),
    _pqs{ds::BinaryMaxHeap<HypernodeID, Gain>(_hg.initialNumNodes()),
         ds::BinaryMaxHeap<HypernodeID, Gain>(_hg.initialNumNodes())},
    _fixtures{std::vector<HypernodeID>(_hg.initialNumNodes()),
              std::vector<HypernodeID>(_hg.initialNumNodes())} {
    _qg = std::make_shared<AdjacencyMatrixQuotientGraph<DFSCycleDetector>>(_hg, _context);
    _gain_manager = std::make_shared<KMinusOneGainManager>(_hg, _context);
    _control_qg_and_gm = true;
  }

  ~AcyclicTwoWayKMinusOneRefiner() override = default;

  AcyclicTwoWayKMinusOneRefiner(const AcyclicTwoWayKMinusOneRefiner&) = delete;
  AcyclicTwoWayKMinusOneRefiner& operator=(const AcyclicTwoWayKMinusOneRefiner&) = delete;

  AcyclicTwoWayKMinusOneRefiner(AcyclicTwoWayKMinusOneRefiner&&) = delete;
  AcyclicTwoWayKMinusOneRefiner& operator=(AcyclicTwoWayKMinusOneRefiner&&) = delete;

  const std::vector<Move>& moves() {
    return _moves;
  }

  void preUncontraction(const HypernodeID representant) override {
    removeFixtures(representant);
    _hns_to_activate.clear(); // populated by removeFixtures()

    if (_control_qg_and_gm) {
      _qg->preUncontraction(representant);
      _gain_manager->preUncontraction(representant);
    }
  }

  void postUncontraction(const HypernodeID representant, const std::vector<HypernodeID>&& partners) override {
    if (_control_qg_and_gm) {
      _qg->postUncontraction(representant, std::forward<const std::vector<HypernodeID>>(partners));
      _gain_manager->postUncontraction(representant, std::forward<const std::vector<HypernodeID>>(partners));
    }

    ASSERT(partners.size() == 1, "Currently only supports pair contractions.");
    const HypernodeID partner = partners.front();

    addFixtures(representant);
    addFixtures(partner);
    initializeFixtures(representant);
    initializeFixtures(partner);

    _hns_to_deactivate.clear(); // populated by addFixtures()
  }

  void printSummary() const override {
    LOG << "[AcyclicTwoWayKMinusOneRefiner] Number of moves:" << _num_moves;
    LOG << "[AcyclicTwoWayKMinusOneRefiner] Number of moves in the last iteration:" << _num_moves_in_last_iteration;
    LOG << "[AcyclicTwoWayKMinusOneRefiner] Positive gain moves:" << _num_positive_gain_moves;
    LOG << "[AcyclicTwoWayKMinusOneRefiner] Zero gain moves:" << _num_zero_gain_moves;
    LOG << "[AcyclicTwoWayKMinusOneRefiner] Negative gain moves:" << _num_negative_gain_moves;
    LOG << "[AcyclicTwoWayKMinusOneRefiner] Improved imbalance by:" << _improved_imbalance;
    LOG << "[AcyclicTwoWayKMinusOneRefiner] Improved KM1 by:" << _improved_km1;
    LOG << "[AcyclicTwoWayKMinusOneRefiner] Init time:" << _time_init;
    LOG << "[AcyclicTwoWayKMinusOneRefiner] Refine time:" << _time_total_refinement;
    LOG << "[AcyclicTwoWayKMinusOneRefiner] Gains time:" << _time_update_gains;
    LOG << "[AcyclicTwoWayKMinusOneRefiner] Rollback time:" << _time_rollback;
    LOG << "[AcyclicTwoWayKMinusOneRefiner] Fixture init+update time:" << _time_init_and_update_fixtures;
    LOG << "[AcyclicTwoWayKMinusOneRefiner] QG update time:" << _time_update_qg;
    LOG << "[AcyclicTwoWayKMinusOneRefiner] HN activation time:" << _time_hn_activate;
    LOG << "[AcyclicTwoWayKMinusOneRefiner] HN deactivation time:" << _time_hn_deactivate;
  }

 private:
  void initializeImpl(const HyperedgeWeight max_gain) final {
    _timer.start();
    _is_initialized = true;

    if (_gain_manager) {
      _qg->rebuild();
      _gain_manager->initialize();
    }
    _gain_manager->resetDelta();
    _hg.resetHypernodeState();

    _num_moves = 0;
    _num_moves_in_last_iteration = 0;
    _num_positive_gain_moves = 0;
    _num_negative_gain_moves = 0;
    _num_zero_gain_moves = 0;
    _improved_imbalance = 0.0;
    _improved_km1 = 0;
    _time_total_refinement = 0.0;
    _time_update_gains = 0.0;
    _time_rollback = 0.0;
    _time_init_and_update_fixtures = 0.0;
    _time_update_qg = 0.0;
    _time_hn_activate = 0.0;
    _time_hn_deactivate = 0.0;

    _fixtures[PREDECESSORS].clear();
    _fixtures[PREDECESSORS].resize(_hg.initialNumNodes());
    _fixtures[SUCCESSORS].clear();
    _fixtures[SUCCESSORS].resize(_hg.initialNumNodes());
    _timer.start();
    initializeFixtures();
    _time_init_and_update_fixtures += _timer.stop();

    _time_init = _timer.stop(); // this should be = and not +=
  }

  void performMovesAndUpdateCacheImpl(const std::vector<Move>& moves,
                                      std::vector<HypernodeID>& refinement_nodes,
                                      const UncontractionGainChanges& changes) final {
    throw std::runtime_error("not implemented");
  }

  bool refineImpl(std::vector<HypernodeID>& refinement_nodes,
                  const std::array<HypernodeWeight, 2>&,
                  const UncontractionGainChanges&,
                  Metrics& best_metrics) final {
    ASSERT(_hg.k() == 2, "2way refiner, but graph has more than 2 blocks:" << _hg.k());
    ASSERT(_context.local_search.fm.max_number_of_fruitless_moves > 0, "Bad configuration, number of fruitless moves is zero");
    DBG << "Refiner call with" << refinement_nodes;

    _timer.start(); // full refinement timer

    ASSERT(best_metrics.imbalance == metrics::imbalance(_hg, _context));
    ASSERT(best_metrics.km1 == metrics::km1(_hg));
    const double initial_imbalance = best_metrics.imbalance;
    const HyperedgeWeight initial_km1 = best_metrics.km1;
    HyperedgeWeight current_km1 = initial_km1;

    _hg.resetHypernodeState();
    _moves.clear();
    _pqs[0].clear();
    _pqs[1].clear();

    ASSERT(VALIDATE_FIXTURES_STATE());
    ASSERT(_qg->isAcyclic());
    ASSERT(AdjacencyMatrixQuotientGraph<DFSCycleDetector>(_hg, _context).isAcyclic());
    ASSERT(!_gain_manager->hasDelta());

    // activate refinement nodes
    for (const HypernodeID& hn : refinement_nodes) {
      if (!_hg.active(hn) && isMovable(hn) && _hg.isBorderNode(hn)) {
        activate(hn);
      }

      for (const HyperedgeID& he : _hg.incidentEdges(hn)) {
        for (const HypernodeID& pin : _hg.pins(he)) {
          if (!_hg.active(pin) && _hg.isBorderNode(pin) && isMovable(pin)) {
            activate(pin);
          }
        }
      }
    }
    DBG << "Initial PQ sizes:" << _pqs[0].size() << "--" << _pqs[1].size();

    // main refinement loop
    int moves_since_improvement = 0;
    int min_cut_index = -1;

    while (moves_since_improvement < 350 // TODO
           && (!_pqs[0].empty() || !_pqs[1].empty())) {
      ASSERT(VALIDATE_PQS_STATE());
      ASSERT(VALIDATE_FIXTURES_STATE());

      // choose from and to part
      const PartitionID from_part = selectPQ();
      ASSERT(from_part != Hypergraph::kInvalidPartition); // already check for that in loop condition
      ASSERT(!_pqs[from_part].empty());
      const PartitionID to_part = otherPartition(from_part);

      // pull hypernode from selected PQ
      const HypernodeID max_gain_hn = _pqs[from_part].top();
      const Gain max_gain = _pqs[from_part].topKey();
      _pqs[from_part].pop();

      ASSERT(_hg.partID(max_gain_hn) == from_part);
      ASSERT(isMovable(max_gain_hn));
      ASSERT(max_gain == _gain_manager->gain(max_gain_hn, to_part));
      ASSERT(max_gain == gainInducedByHypergraph(max_gain_hn, to_part));
      ASSERT(_hg.active(max_gain_hn));
      ASSERT(!_hg.marked(max_gain_hn));

      _hg.mark(max_gain_hn); // exclude this hypernode from this round

      // check whether we can move without violating the imbalance constraint
      if (_hg.partWeight(to_part) + _hg.nodeWeight(max_gain_hn) <= _context.partition.max_part_weights[to_part]) {
        // gather some statistics
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

        // perform actual hypernode movement
        current_km1 -= max_gain;
        if (_hg.nodeIsEnabled(hn_to_debug)) DBG << "Before move:" << isTailPartition(0) << isHeadPartition(0) << isTailPartition(1) << isHeadPartition(1) << isMovable(hn_to_debug) << fixtures(hn_to_debug) << _fixtures[0][hn_to_debug] << _fixtures[1][hn_to_debug];
        DBG << "Move HN" << max_gain_hn << ":" << from_part << "-->" << to_part << "for" << max_gain << ","
            << initial_km1 << "-->" << best_metrics.km1 << "-->" << current_km1;
        move(max_gain_hn, from_part, to_part);
        if (_hg.nodeIsEnabled(hn_to_debug)) DBG << "After move:" << isTailPartition(0) << isHeadPartition(0) << isTailPartition(1) << isHeadPartition(1) << isMovable(hn_to_debug) << fixtures(hn_to_debug) << _fixtures[0][hn_to_debug] << _fixtures[1][hn_to_debug];

        ++moves_since_improvement;
        ASSERT(current_km1 == metrics::km1(_hg));

        // accept new partition if it improves the cut
        bool accept = current_km1 < best_metrics.km1;
        if (!accept && current_km1 == best_metrics.km1) {
          accept = Randomize::instance().flipCoin();
        }
        if (accept) {
          DBG << "Accept improvement:" << best_metrics.km1 << "-->" << current_km1;
          best_metrics.km1 = current_km1;
          _gain_manager->resetDelta();
          moves_since_improvement = 0;
          min_cut_index = _moves.size();
        }

        _moves.emplace_back(max_gain_hn, from_part, to_part);
      } else { // move infeasible due to imbalance constraint, but we can still activate its neighbors
        for (const HyperedgeID& he : _hg.incidentEdges(max_gain_hn)) {
          for (const HypernodeID& hn : _hg.pins(he)) {
            if (!_hg.marked(hn) && !_hg.active(hn) && _hg.isBorderNode(hn) && isMovable(hn)) {
              activate(hn);
            }
          }
        }
      }
    }

    ASSERT(VALIDATE_FIXTURES_STATE());

    // rollback to the best last accepted state
    _timer.start();
    int last_index = _moves.size() - 1;
    while (last_index != min_cut_index) {
      const HypernodeID hn = _moves[last_index].hn;
      const PartitionID from_part = _moves[last_index].to;
      const PartitionID to_part = _moves[last_index].from;
      _moves.pop_back();
      --last_index;

      // revert hypernode movement
      DBG << "Rollback" << hn << ":" << from_part << "-->" << to_part;
      removeFixtures(hn);
      const bool success = _qg->testAndUpdateBeforeMovement(hn, from_part, to_part);
      ASSERT(success);
      _hg.changeNodePart(hn, from_part, to_part);
      addFixtures(hn);
      initializeFixtures(hn);

      // statistics
      --_num_moves;
      if (_hg.currentNumNodes() == _hg.initialNumNodes()) {
        --_num_moves_in_last_iteration;
      }
    }
    _hns_to_activate.clear();
    _hns_to_deactivate.clear();
    _gain_manager->rollbackDelta();
    _time_rollback += _timer.stop();

    ASSERT(_qg->isAcyclic());
    ASSERT(AdjacencyMatrixQuotientGraph<DFSCycleDetector>(_hg, _context).isAcyclic());
    ASSERT(!_gain_manager->hasDelta());
    ASSERT(best_metrics.km1 == metrics::km1(_hg));

    best_metrics.imbalance = metrics::imbalance(_hg, _context);
    _improved_imbalance += initial_imbalance - best_metrics.imbalance;
    _improved_km1 += initial_km1 - best_metrics.km1;

    ASSERT(_timer.running());
    _time_total_refinement += _timer.stop();
    ASSERT(!_timer.running());

    return CutDecreasedOrInfeasibleImbalanceDecreased::improvementFound(best_metrics.km1, initial_km1,
                                                                        best_metrics.imbalance, initial_imbalance,
                                                                        _context.partition.epsilon);
  }

  PartitionID selectPQ() const {
    // handle empty PQs
    if (_pqs[0].empty() && _pqs[1].empty()) {
      return Hypergraph::kInvalidPartition;
    }
    if (_pqs[0].empty()) {
      return 1;
    }
    if (_pqs[1].empty()) {
      return 0;
    }

    // both PQs viable, select overloaded block if any
    if (isOverloaded(0)) {
      return 0;
    }
    if (isOverloaded(1)) {
      return 1;
    }

    // TODO this seems to work *A LOT* better than selecting the PQ based on topKey()
    if (_hg.partWeight(0) > _hg.partWeight(1)) {
      return 0;
    } else {
      return 1;
    }

    // both PQs viable and partition is balanced, select based on gain
    if (_pqs[0].topKey() > _pqs[1].topKey()) {
      return 0;
    }
    if (_pqs[1].topKey() > _pqs[0].topKey()) {
      return 1;
    }
    ASSERT(_pqs[0].topKey() == _pqs[1].topKey());
    return Randomize::instance().flipCoin() ? 1 : 0; // both PQs offer the same gain so flip a coin to decide
  }

  bool isOverloaded(const PartitionID part) const {
    return _hg.partWeight(part) > _context.partition.max_part_weights[part];
  }

  void move(const HypernodeID hn, const PartitionID from_part, const PartitionID to_part) {
    ASSERT(_hg.marked(hn), "Mark HN before calling move()");
    ASSERT(_hg.partID(hn) == from_part);

    ASSERT(_hns_to_activate.empty());
    ASSERT(_hns_to_deactivate.empty());

#ifdef KAHYPAR_USE_ASSERTIONS
    const auto km1_before_move = metrics::km1(_hg);
    const auto expected_gain = _gain_manager->gain(hn, to_part);
#endif // KAHYPAR_USE_ASSERTIONS

    _timer.start();
    removeFixtures(hn);
    _time_init_and_update_fixtures += _timer.stop();

    _timer.start();
    const bool success = _qg->testAndUpdateBeforeMovement(hn, from_part, to_part);
    ASSERT(success, V(hn) << V(from_part) << V(to_part));
    _time_update_qg += _timer.stop();

    _hg.changeNodePart(hn, from_part, to_part);

    _timer.start();
    _gain_manager->updateAfterMovement(hn, from_part, to_part,
                                       [&](const HypernodeID& hn, const PartitionID& part, const Gain& gain) {
                                         if (!_hg.marked(hn) && !_hg.active(hn) && isMovable(hn)) {
                                           _hns_to_activate.push_back(hn);
                                           return;
                                         }
                                       },
                                       [&](const HypernodeID& hn, const PartitionID& part, const Gain& delta) {
                                         if (_hg.marked(hn) || !_hg.active(hn) ||
                                             part == Hypergraph::kInvalidPartition) {
                                           return;
                                         }

                                         const auto source_part = _hg.partID(hn);
                                         ASSERT(_pqs[source_part].contains(hn));
                                         _pqs[source_part].updateKeyBy(hn, delta);
                                       },
                                       [&](const HypernodeID& hn, const PartitionID& part, const Gain&) {
                                         if (!_hg.marked(hn) && _hg.active(hn)) {
                                           _hns_to_deactivate.push_back(hn);
                                         }
                                       });
    _time_update_gains += _timer.stop();

    // activate new movable nodes adjacent to `hn`
#ifdef KAHYPAR_USE_ASSERTIONS
    for (const HypernodeID& hn_to_activate : _hns_to_activate) {
      ASSERT(!_hg.marked(hn_to_activate));
      ASSERT(!_hg.active(hn_to_activate));
    }
#endif // KAHYPAR_USE_ASSERTIONS
    for (const HypernodeID& hn_to_activate : _hns_to_activate) {
      if (_hg.active(hn_to_activate)) { // hn might be contained multiple times in _hns_to_activate
        continue;
      }
      activate(hn_to_activate);
    }
    _hns_to_activate.clear(); // populated by removeFixtures() and during gain update

    _timer.start();
    initializeFixtures(hn);
    addFixtures(hn);
    _time_init_and_update_fixtures += _timer.stop();

    // deactivate nodes that are no longer movable
#ifdef KAHYPAR_USE_ASSERTIONS
    for (const HypernodeID& hn_to_deactivate : _hns_to_deactivate) {
      ASSERT(!_hg.marked(hn_to_deactivate));
      ASSERT(_hg.active(hn_to_deactivate));
    }
#endif // KAHYPAR_USE_ASSERTIONS
    for (const HypernodeID& hn_to_deactivate : _hns_to_deactivate) {
      DBGC(hn_to_deactivate == hn_to_debug) << "Deactivate" << hn_to_deactivate << _hg.active(hn_to_deactivate);
      if (!_hg.active(hn_to_deactivate)) { // hn might be contained multiple times in _hns_to_deactivate
        continue;
      }
      deactivate(hn_to_deactivate);
    }
    _hns_to_deactivate.clear(); // populated by addFixtures() and during gain update

#ifdef KAHYPAR_USE_ASSERTIONS
    const auto km1_after_move = metrics::km1(_hg);
    ASSERT(km1_before_move - expected_gain == km1_after_move,
           "Before:" << V(km1_before_move) << "Gain:" << V(expected_gain) << "After:" << V(km1_after_move));
#endif // KAHYPAR_USE_ASSERTIONS
  }

  void activate(const HypernodeID hn) {
    ASSERT(isMovable(hn));
    ASSERT(_hg.isBorderNode(hn));
    ASSERT(!_hg.active(hn));
    ASSERT(!_hg.marked(hn));

    _timer.start();
    const PartitionID source = _hg.partID(hn);
    const PartitionID target = otherPartition(source);
    const Gain gain = _gain_manager->gain(hn, target);
    _pqs[source].push(hn, gain);
    _hg.activate(hn);
    _time_hn_activate += _timer.stop();

    DBG << "Activated" << hn << "(" << _hg.partID(hn) << ") with gain" << gain;
  }

  void deactivate(const HypernodeID hn) {
    ASSERT(!_hg.marked(hn));
    ASSERT(_hg.active(hn));

    _timer.start();
    const auto part = _hg.partID(hn);
    ASSERT(_pqs[part].contains(hn));
    _pqs[part].remove(hn);
    _hg.deactivate(hn);
    _time_hn_deactivate += _timer.stop();

    DBG << "Deactivate" << hn << "(" << _hg.partID(hn) << ")";
  }

  void removeFixtures(const HypernodeID hn) {
    for (const HyperedgeID& he : _hg.incidentTailEdges(hn)) {
      for (const HypernodeID& head : _hg.heads(he)) {
        if (_hg.partID(head) == _hg.partID(hn)) {
          ASSERT(_fixtures[PREDECESSORS][head] > 0);
          --_fixtures[PREDECESSORS][head];

          if (!_hg.active(head) && isMovable(head)) {
            _hns_to_activate.push_back(head);
          }
        }
      }
    }

    for (const HyperedgeID& he : _hg.incidentHeadEdges(hn)) {
      for (const HypernodeID& tail : _hg.tails(he)) {
        if (_hg.partID(tail) == _hg.partID(hn)) {
          ASSERT(_fixtures[SUCCESSORS][tail] > 0);
          --_fixtures[SUCCESSORS][tail];

          if (!_hg.active(tail) && isMovable(tail)) {
            _hns_to_activate.push_back(tail);
          }
        }
      }
    }
  }

  void addFixtures(const HypernodeID hn) {
    for (const HyperedgeID& he : _hg.incidentTailEdges(hn)) {
      for (const HypernodeID& head : _hg.heads(he)) {
        if (_hg.partID(head) == _hg.partID(hn)) {
          DBGC(head == hn_to_debug) << "Increment _fixtures[" << 1 << "][" << head << "]";
          ++_fixtures[PREDECESSORS][head];

          if (_hg.active(head) && !isMovable(head)) {
            _hns_to_deactivate.push_back(head);
          }
        }
      }
    }

    for (const HyperedgeID& he : _hg.incidentHeadEdges(hn)) {
      for (const HypernodeID& tail : _hg.tails(he)) {
        if (_hg.partID(tail) == _hg.partID(hn)) {
          DBGC(tail == hn_to_debug) << "Increment _fixtures[" << 0 << "][" << tail << "]";
          ++_fixtures[SUCCESSORS][tail];

          if (_hg.active(tail) && !isMovable(tail)) {
            _hns_to_deactivate.push_back(tail);
          }
        }
      }
    }
  }

  void initializeFixtures() {
    for (const HypernodeID& hn : _hg.nodes()) {
      initializeFixtures(hn);
    }
  }

  void initializeFixtures(const HypernodeID hn) {
    _fixtures[SUCCESSORS][hn] = 0;
    _fixtures[PREDECESSORS][hn] = 0;

    for (const HyperedgeID& he : _hg.incidentTailEdges(hn)) {
      for (const HypernodeID& head : _hg.heads(he)) {
        if (_hg.partID(head) == _hg.partID(hn)) {
          DBGC(hn == hn_to_debug) << "Increment _fixtures[" << 0 << "][" << hn << "]";
          ++_fixtures[SUCCESSORS][hn];

          if (_hg.active(hn) && !isMovable(hn)) {
            _hns_to_deactivate.push_back(hn);
          }
        }
      }
    }

    for (const HyperedgeID& he : _hg.incidentHeadEdges(hn)) {
      for (const HypernodeID& tail : _hg.tails(he)) {
        if (_hg.partID(tail) == _hg.partID(hn)) {
          DBGC(hn == hn_to_debug) << "Increment _fixtures[" << 1 << "][" << hn << "]";
          ++_fixtures[PREDECESSORS][hn];

          if (_hg.active(hn) && !isMovable(hn)) {
            _hns_to_deactivate.push_back(hn);
          }
        }
      }
    }
  }

  bool isMovable(const HypernodeID hn) const {
    return fixtures(hn) == 0;
  }

  HypernodeID fixtures(const HypernodeID hn) const {
    const PartitionID source = _hg.partID(hn);
    return isTailPartition(source)
           ? _fixtures[SUCCESSORS][hn]
           : _fixtures[PREDECESSORS][hn];
  }

  bool isTailPartition(const PartitionID part) const {
    if (_qg->connected(part, 1 - part)) {
      return true;
    } else if (!_qg->connected(1 - part, part)) {
      // catch edge case where there are no edges in the quotient graph: is this case,
      // make part 0 the tail and part 1 the head
      return part == 0;
    }
    return false;
  }

  bool isHeadPartition(const PartitionID part) const {
    return !isTailPartition(part);
  }

  static PartitionID otherPartition(const PartitionID part) {
    return 1 - part;
  }

#ifdef KAHYPAR_USE_ASSERTIONS

  Gain gainInducedByHypergraph(const HypernodeID hn, const PartitionID target_part) const {
    Gain gain = 0;
    for (const HyperedgeID& he : _hg.incidentEdges(hn)) {
      //ASSERT(_hg.edgeSize(he) > 1, V(he));
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

  bool VALIDATE_FIXTURES_STATE() const {
    for (const HypernodeID& hn : _hg.nodes()) {
      HypernodeID num_predecessors = 0;
      HypernodeID num_successors = 0;
      for (const HyperedgeID& he : _hg.incidentTailEdges(hn)) {
        for (const HypernodeID& head : _hg.heads(he)) {
          if (_hg.partID(head) == _hg.partID(hn)) {
            ++num_successors;
          }
        }
      }
      for (const HyperedgeID& he : _hg.incidentHeadEdges(hn)) {
        for (const HypernodeID& tail : _hg.tails(he)) {
          if (_hg.partID(tail) == _hg.partID(hn)) {
            ++num_predecessors;
          }
        }
      }

      ASSERT(_fixtures[PREDECESSORS][hn] == num_predecessors);
      ASSERT(_fixtures[SUCCESSORS][hn] == num_successors);
    }

    return true;
  }

  bool VALIDATE_PQS_STATE() const {
    for (const HypernodeID& hn : _hg.nodes()) {
      if (_hg.marked(hn)) {
        ASSERT(!_pqs[0].contains(hn));
        ASSERT(!_pqs[1].contains(hn));
        continue;
      }

      if (_hg.active(hn)) {
        const auto source_part = _hg.partID(hn);
        const auto target_part = otherPartition(source_part);
        //ASSERT(_gain_manager.isAdjacentTo(hn, target_part));
        ASSERT(_pqs[source_part].contains(hn));
        ASSERT(!_pqs[target_part].contains(hn));
        ASSERT(_pqs[source_part].getKey(hn) == _gain_manager->gain(hn, target_part));
        ASSERT(_pqs[source_part].getKey(hn) == gainInducedByHypergraph(hn, target_part));
        ASSERT(isMovable(hn));
        continue;
      }

      // not active, not marked
      ASSERT(!_pqs[0].contains(hn));
      ASSERT(!_pqs[1].contains(hn));
    }
    return true;
  }

#endif

  Hypergraph& _hg;
  const Context& _context;
  std::shared_ptr<AdjacencyMatrixQuotientGraph<DFSCycleDetector>> _qg;
  std::shared_ptr<KMinusOneGainManager> _gain_manager;
  std::array<ds::BinaryMaxHeap<HypernodeID, Gain>, 2> _pqs;
  std::array<std::vector<HypernodeID>, 2> _fixtures;
  std::vector<Move> _moves;
  std::vector<HypernodeID> _hns_to_activate;
  std::vector<HypernodeID> _hns_to_deactivate;

  std::size_t _num_moves{0};
  std::size_t _num_positive_gain_moves{0};
  std::size_t _num_negative_gain_moves{0};
  std::size_t _num_zero_gain_moves{0};
  std::size_t _num_moves_in_last_iteration{0};
  HyperedgeWeight _improved_km1{0};
  double _improved_imbalance{0.0};
  double _time_total_refinement{0.0};
  double _time_rollback{0.0};
  double _time_update_gains{0.0};
  double _time_init_and_update_fixtures{0.0};
  double _time_update_qg{0.0};
  double _time_hn_activate{0.0};
  double _time_hn_deactivate{0.0};
  double _time_init{0.0};
  bool _control_qg_and_gm{false};

  HTimer _timer;
};
} // namespace kahypar