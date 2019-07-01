#pragma once

#include "kahypar/datastructure/hypergraph.h"
#include "kahypar/partition/refinement/i_refiner.h"
#include "kahypar/partition/refinement/acyclic_soft_rebalance_refiner.h"
#include "kahypar/partition/refinement/acyclic_kway_fm_km1_refiner.h"
#include "kahypar/partition/refinement/acyclic_hard_rebalance_refiner.h"
#include "kahypar/partition/refinement/policies/fm_stop_policy.h"
#include "kahypar/partition/context.h"

namespace kahypar {
class AcyclicKMinusOneRefiner final : public IRefiner {
 private:
  constexpr static bool debug = false;

 public:
  AcyclicKMinusOneRefiner(Hypergraph& hypergraph, const Context& context) :
    _hg(hypergraph),
    _original_context(context),
    _local_search(hypergraph, context),
    _soft_rebalance(hypergraph, context, _local_search.qg()),
    _hard_rebalance(hypergraph, context, _local_search.qg()) {}

  ~AcyclicKMinusOneRefiner() override = default;

  AcyclicKMinusOneRefiner(const AcyclicKMinusOneRefiner&) = delete;
  AcyclicKMinusOneRefiner& operator=(const AcyclicKMinusOneRefiner&) = delete;

  AcyclicKMinusOneRefiner(AcyclicKMinusOneRefiner&&) = delete;
  AcyclicKMinusOneRefiner& operator=(AcyclicKMinusOneRefiner&&) = delete;

  void preUncontraction(const HypernodeID u) override {
    _local_search.qg().preUncontraction(u);
    ASSERT(dag::isAcyclic(_hg), V(u));
  }

  void postUncontraction(const HypernodeID u, const HypernodeID v) override {
    _local_search.qg().postUncontraction(u, v);
    ASSERT(dag::isAcyclic(_hg), V(u) << V(v));
  }

  void printFinalInfo() override {
    LOG << "Moves by _local_search:" << _local_search.numMoves();
    LOG << "Moves by _soft_rebalance:" << _soft_rebalance.numMoves();
    LOG << "Imbalance improved by _soft_rebalance:" << _soft_rebalance.improvedImbalanceBy();
    LOG << "Rejected due to balance:" << _soft_rebalance.numRejectedDueToBalance();
    LOG << "Rejected due to KM1:" << _soft_rebalance.numRejectedDueToKM1();
    LOG << "Imbalance improved by _hard_rebalance:" << _hard_rebalance.improvedImbalance();
    LOG << "Positive Gains in _soft_rebalance:" << _soft_rebalance.numPosGains();
    LOG << "Zero Gains in _soft_rebalance:" << _soft_rebalance.numZeroGains();
    LOG << "Negative Gains in _soft_rebalance:" << _soft_rebalance.numNegGains();
    LOG << "Positive Gains in _local_search:" << _local_search.numPosGains();
    LOG << "Zero Gains in _local_search:" << _local_search.numZeroGains();
    LOG << "Negative Gains in _local_search:" << _local_search.numNegGains();
  }

 private:
  void initializeImpl(const HyperedgeWeight max_gain) final {
    if (_original_context.enable_soft_rebalance) {
      _soft_rebalance.initialize(max_gain);
    }
    _local_search.initialize(max_gain);
    _hard_rebalance.initialize(max_gain);
    _least_num_nodes = _hg.currentNumNodes();
    _is_initialized = true;
  }

  void performMovesAndUpdateCacheImpl(const std::vector<Move>& moves,
                                      std::vector<HypernodeID>& refinement_nodes,
                                      const UncontractionGainChanges& changes) final {
    if (_original_context.enable_soft_rebalance) {
      _soft_rebalance.performMovesAndUpdateCache(moves, refinement_nodes, changes);
    }
    _local_search.performMovesAndUpdateCache(moves, refinement_nodes, changes);
    _hard_rebalance.performMovesAndUpdateCache(moves, refinement_nodes, changes);
  }

  bool refineImpl(std::vector<HypernodeID>& refinement_nodes,
                  const std::array<HypernodeWeight, 2>& max_allowed_part_weights,
                  const UncontractionGainChanges& uncontraction_changes,
                  Metrics& best_metrics) final {
    double progress = 1.0 - (_hg.initialNumNodes() - _hg.currentNumNodes()) /
                            static_cast<double>(_hg.initialNumNodes() - _least_num_nodes);
    double delta_epsilon = _original_context.partition.epsilon - _original_context.partition.final_epsilon;
    double epsilon = _original_context.partition.epsilon - delta_epsilon * progress;
    const double current_imbalance = best_metrics.imbalance;
    DBG << V(epsilon) << V(current_imbalance);

    if (_original_context.enable_soft_rebalance) {
      _soft_rebalance.changeEpsilon(epsilon);
    }
    _hard_rebalance.changeEpsilon(epsilon);
    _local_search.changeEpsilon(epsilon);

    _hard_rebalance.enableRefinementNodes(refinement_nodes);
    if (_original_context.enable_soft_rebalance) {
      _soft_rebalance.enableRefinementNodes(refinement_nodes);
    }

    if (_original_context.enable_soft_rebalance && current_imbalance > epsilon) {
      _soft_rebalance.refine(refinement_nodes, max_allowed_part_weights, uncontraction_changes, best_metrics);
      const auto& moves = mergeMultiMoves(_soft_rebalance.moves());
      _hard_rebalance.performMovesAndUpdateCache(moves, refinement_nodes, uncontraction_changes);
      _local_search.performMovesAndUpdateCache(moves, refinement_nodes, uncontraction_changes);
      DBG << "_soft_rebalance" << current_imbalance << "-->" << best_metrics.imbalance << V(epsilon);
    }

    if (current_imbalance > epsilon) {
      _hard_rebalance.refine(refinement_nodes, max_allowed_part_weights, uncontraction_changes, best_metrics);
      const auto& moves = mergeMultiMoves(_hard_rebalance.moves());
      if (_original_context.enable_soft_rebalance) {
        _soft_rebalance.performMovesAndUpdateCache(moves, refinement_nodes, uncontraction_changes);
      }
      _local_search.performMovesAndUpdateCache(moves, refinement_nodes, uncontraction_changes);
      DBG << "_hard_rebalance:" << current_imbalance << "-->" << best_metrics.imbalance << V(epsilon);
    }

    bool improvement = _local_search.refine(refinement_nodes, max_allowed_part_weights, uncontraction_changes, best_metrics);
    const auto& moves = _local_search.moves();
    if (_original_context.enable_soft_rebalance) {
      _soft_rebalance.performMovesAndUpdateCache(moves, refinement_nodes, uncontraction_changes);
    }
    _hard_rebalance.performMovesAndUpdateCache(moves, refinement_nodes, uncontraction_changes);

    return improvement;
  }

  std::vector<Move> mergeMultiMoves(const std::vector<Move>& moves) {
    std::vector<std::size_t> last_move(_hg.initialNumNodes());
    std::vector<std::size_t> first_move(_hg.initialNumNodes());
    for (std::size_t i = moves.size(); i > 0; --i) {
      const Move& move = moves[i - 1];
      ASSERT(move.hn < _hg.initialNumNodes());
      first_move[move.hn] = i - 1;
    }
    for (std::size_t i = 0; i < moves.size(); ++i) {
      last_move[moves[i].hn] = i;
    }
    std::vector<Move> merged_moves;
    for (std::size_t i = 0; i < moves.size(); ++i) {
      const Move& move = moves[i];
      if (i == first_move[move.hn] && move.from != moves[last_move[move.hn]].to) {
        merged_moves.emplace_back(move.hn, move.from, moves[last_move[move.hn]].to);
      }
    }

    for (const Move& move : merged_moves) {
      ASSERT(_hg.partID(move.hn) == move.to);
    }

    return merged_moves;
  }

  /*********************************************************************************
   * MEMBER FIELDS
   *********************************************************************************/
  Hypergraph& _hg;
  const Context& _original_context;

  AcyclicKWayKMinusOneRefiner<NumberOfFruitlessMovesStopsSearch> _local_search;
  AcyclicSoftRebalanceRefiner _soft_rebalance;
  AcyclicHardRebalanceRefiner _hard_rebalance;

  HypernodeID _least_num_nodes;
};

} // namespace kahypar