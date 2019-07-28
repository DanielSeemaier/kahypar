#pragma once

#include "kahypar/datastructure/hypergraph.h"
#include "kahypar/partition/refinement/i_refiner.h"
#include "kahypar/partition/refinement/acyclic_local_search_refiner.h"
#include "kahypar/partition/refinement/acyclic_hard_rebalance_refiner.h"
#include "kahypar/partition/refinement/acyclic_soft_rebalance_refiner.h"
#include "kahypar/partition/refinement/policies/fm_stop_policy.h"
#include "kahypar/partition/context.h"

#include "kahypar/partition/refinement/km1_gain_manager.h"

namespace kahypar {
class AcyclicKMinusOneRefiner final : public IRefiner {
 private:
  constexpr static bool debug = false;

 public:
  AcyclicKMinusOneRefiner(Hypergraph& hypergraph, const Context& context) :
    _hg(hypergraph),
    _original_context(context),
    _context(context),
    _qg(hypergraph, context),
    _gain_manager(hypergraph, context),
    _local_search_refiner(hypergraph, _context, _qg, _gain_manager),
    _hard_rebalance_refiner(hypergraph, _context, _qg, _gain_manager),
    _soft_rebalance_refiner(hypergraph, _context, _qg, _gain_manager) {}

  ~AcyclicKMinusOneRefiner() override = default;

  AcyclicKMinusOneRefiner(const AcyclicKMinusOneRefiner&) = delete;
  AcyclicKMinusOneRefiner& operator=(const AcyclicKMinusOneRefiner&) = delete;

  AcyclicKMinusOneRefiner(AcyclicKMinusOneRefiner&&) = delete;
  AcyclicKMinusOneRefiner& operator=(AcyclicKMinusOneRefiner&&) = delete;

  void preUncontraction(const HypernodeID representant) override {
    _local_search_refiner.preUncontraction(representant);
    _soft_rebalance_refiner.preUncontraction(representant);
    _hard_rebalance_refiner.preUncontraction(representant);
    _gain_manager.preUncontraction(representant);
    _qg.preUncontraction(representant);
  }

  void postUncontraction(const HypernodeID representant, const std::vector<HypernodeID>&& partners) override {
    _qg.postUncontraction(representant, std::forward<const std::vector<HypernodeID>>(partners));
    _gain_manager.postUncontraction(representant, std::forward<const std::vector<HypernodeID>>(partners));
    _hard_rebalance_refiner.postUncontraction(representant, std::forward<const std::vector<HypernodeID>>(partners));
    _soft_rebalance_refiner.postUncontraction(representant, std::forward<const std::vector<HypernodeID>>(partners));
    _local_search_refiner.postUncontraction(representant, std::forward<const std::vector<HypernodeID>>(partners));
  }

  void printSummary() const override {
    if (_context.enable_soft_rebalance) {
      LOG << "Soft Rebalance Refiner:";
      _soft_rebalance_refiner.printSummary();
    } else {
      LOG << "Soft Rebalance disabled";
    }
    LOG << "Hard Rebalance Refiner:";
    _hard_rebalance_refiner.printSummary();
    LOG << "Local Search Refiner:";
    _local_search_refiner.printSummary();
    LOG << "Running with extra refinement nodes:" << _context.refine_rebalance_moves;
  }

 private:
  void initializeImpl(const HyperedgeWeight max_gain) final {
    _qg.rebuild();
    _gain_manager.initialize();
    _local_search_refiner.initialize(max_gain);
    _hard_rebalance_refiner.initialize(max_gain);
    if (_context.enable_soft_rebalance) {
      _soft_rebalance_refiner.initialize(max_gain);
    }
    _least_num_nodes = _hg.currentNumNodes();
    _is_initialized = true;
  }

  void performMovesAndUpdateCacheImpl(const std::vector<Move>& moves,
                                      std::vector<HypernodeID>& refinement_nodes,
                                      const UncontractionGainChanges& changes) final {
    ASSERT(false, "Currently not implemented.");
  }

  double calculateEpsilon() const {
    if (_least_num_nodes == _hg.initialNumNodes()) { // special case: refine finest graph
      return _original_context.partition.epsilon;
    }

    const double progress = 1.0 - (_hg.initialNumNodes() - _hg.currentNumNodes())
                                  / static_cast<double>(_hg.initialNumNodes() - _least_num_nodes);
    const double delta = _original_context.partition.epsilon - _original_context.partition.final_epsilon;
    return _original_context.partition.epsilon - delta * progress;
  }

  void updateEpsilon() {
    _context.partition.epsilon = calculateEpsilon();
    _context.setupPartWeights(_hg.totalWeight());
  }

  bool refineImpl(std::vector<HypernodeID>& refinement_nodes,
                  const std::array<HypernodeWeight, 2>& max_allowed_part_weights,
                  const UncontractionGainChanges& uncontraction_changes,
                  Metrics& best_metrics) final {
    ASSERT(best_metrics.km1 == metrics::km1(_hg));
    ASSERT(best_metrics.imbalance == metrics::imbalance(_hg, _context));

    std::vector<HypernodeID> more_refinement_nodes = refinement_nodes;

    updateEpsilon();
    const double current_imbalance = best_metrics.imbalance;
    DBG << "Imbalance:" << current_imbalance << "(current) -->" << _context.partition.epsilon << "(goal) with"
        << _hg.currentNumNodes() << "/" << _hg.initialNumNodes() << "nodes";

    if (_context.enable_soft_rebalance && best_metrics.imbalance > _context.partition.epsilon) {
      DBG << "Running soft rebalance because" << best_metrics.imbalance << ">" << _context.partition.epsilon;
      _soft_rebalance_refiner.refine(refinement_nodes, max_allowed_part_weights, uncontraction_changes, best_metrics);
      ASSERT(!_gain_manager.hasDelta());

      const auto& moves = _soft_rebalance_refiner.moves();
      _hard_rebalance_refiner.performMovesAndUpdateCache(moves, refinement_nodes, uncontraction_changes);
      _local_search_refiner.performMovesAndUpdateCache(moves, refinement_nodes, uncontraction_changes);

      if (_context.refine_rebalance_moves) {
        for (const auto& move : moves) {
          if (std::find(more_refinement_nodes.begin(), more_refinement_nodes.end(), move.hn) == more_refinement_nodes.end()) {
            more_refinement_nodes.push_back(move.hn);
          }
        }
      }
      DBG << "-->" << best_metrics.imbalance;
    }

    if (best_metrics.imbalance > _context.partition.epsilon) {
      DBG << "Running hard rebalance because" << best_metrics.imbalance << ">" << _context.partition.epsilon;
      _hard_rebalance_refiner.refine(refinement_nodes, max_allowed_part_weights, uncontraction_changes, best_metrics);
      ASSERT(!_gain_manager.hasDelta());

      const auto& moves = _hard_rebalance_refiner.moves();
      if (_context.enable_soft_rebalance) {
        _soft_rebalance_refiner.performMovesAndUpdateCache(moves, refinement_nodes, uncontraction_changes);
      }
      _local_search_refiner.performMovesAndUpdateCache(moves, refinement_nodes, uncontraction_changes);

      if (_context.refine_rebalance_moves) {
        for (const auto& move : moves) {
          if (std::find(more_refinement_nodes.begin(), more_refinement_nodes.end(), move.hn) == more_refinement_nodes.end()) {
            more_refinement_nodes.push_back(move.hn);
          }
        }
      }
      DBG << "-->" << best_metrics.imbalance;
    }

    const bool stop = _local_search_refiner.refine(more_refinement_nodes, max_allowed_part_weights, uncontraction_changes, best_metrics);
    ASSERT(!_gain_manager.hasDelta());

    const auto& moves = _local_search_refiner.moves();
    if (_context.enable_soft_rebalance) {
      _soft_rebalance_refiner.performMovesAndUpdateCache(moves, refinement_nodes, uncontraction_changes);
    }
    _hard_rebalance_refiner.performMovesAndUpdateCache(moves, refinement_nodes, uncontraction_changes);
    _qg.resetChanged();

    ASSERT(best_metrics.km1 == metrics::km1(_hg));
    ASSERT(best_metrics.imbalance == metrics::imbalance(_hg, _context));
    return stop;
  }

  Hypergraph& _hg;
  const Context& _original_context;
  Context _context;
  HypernodeID _least_num_nodes{0};
  AdjacencyMatrixQuotientGraph<DFSCycleDetector> _qg;
  KMinusOneGainManager _gain_manager;

  AcyclicLocalSearchRefiner<NumberOfFruitlessMovesStopsSearch> _local_search_refiner;
  AcyclicHardRebalanceRefiner _hard_rebalance_refiner;
  AcyclicSoftRebalanceRefiner _soft_rebalance_refiner;
};

} // namespace kahypar