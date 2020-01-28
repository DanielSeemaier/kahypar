#pragma once

#include "kahypar/datastructure/hypergraph.h"
#include "kahypar/partition/refinement/i_refiner.h"
#include "kahypar/partition/refinement/acyclic_local_search_refiner.h"
#include "kahypar/partition/refinement/acyclic_local_search_repeated_refiner.h"
#include "kahypar/partition/refinement/acyclic_hard_rebalance_refiner.h"
#include "kahypar/partition/refinement/acyclic_soft_rebalance_refiner.h"
#include "kahypar/partition/refinement/acyclic_kway_am_fm_refiner.h"
#include "kahypar/partition/refinement/policies/fm_stop_policy.h"
#include "kahypar/partition/context.h"

#include "kahypar/partition/refinement/km1_gain_manager.h"

namespace kahypar {
template<class StoppingPolicy = Mandatory,
    class FMImprovementPolicy = CutDecreasedOrInfeasibleImbalanceDecreased>
class AcyclicKMinusOneRefiner final : public IRefiner {
private:
  constexpr static bool debug = false;

public:
  AcyclicKMinusOneRefiner(Hypergraph &hypergraph, const Context &context) :
      _hg(hypergraph),
      _original_context(context),
      _context(context),
      _qg(std::make_shared<AdjacencyMatrixQuotientGraph<DFSCycleDetector>>(_hg, _context)),
      _gain_manager(std::make_shared<KMinusOneGainManager>(_hg, _context)),
      _local_search_refiner(hypergraph, _context, *_qg, *_gain_manager),
      _local_search_refiner_am(hypergraph, _context, *_qg, *_gain_manager),
      _hard_rebalance_refiner(hypergraph, _context, *_qg, *_gain_manager),
      _soft_rebalance_refiner(hypergraph, _context, *_qg, *_gain_manager) {
  }

  ~AcyclicKMinusOneRefiner() override = default;

  AcyclicKMinusOneRefiner(const AcyclicKMinusOneRefiner &) = delete;

  AcyclicKMinusOneRefiner &operator=(const AcyclicKMinusOneRefiner &) = delete;

  AcyclicKMinusOneRefiner(AcyclicKMinusOneRefiner &&) = delete;

  AcyclicKMinusOneRefiner &operator=(AcyclicKMinusOneRefiner &&) = delete;

  void preUncontraction(const HypernodeID representant) override {
    if (_context.only_do_advanced_moves) {
      _local_search_refiner_am.preUncontraction(representant);
    } else {
      _local_search_refiner.preUncontraction(representant);
    }
    if (_context.enable_soft_rebalance) {
      _soft_rebalance_refiner.preUncontraction(representant);
    }
    if (_context.enable_hard_rebalance) {
      _hard_rebalance_refiner.preUncontraction(representant);
    }
    _gain_manager->preUncontraction(representant);
    _qg->preUncontraction(representant);
  }

  void postUncontraction(const HypernodeID representant, const std::vector<HypernodeID> &&partners) override {
    _qg->postUncontraction(representant, std::forward<const std::vector<HypernodeID>>(partners));
    _gain_manager->postUncontraction(representant, std::forward<const std::vector<HypernodeID>>(partners));
    if (_context.enable_hard_rebalance) {
      _hard_rebalance_refiner.postUncontraction(representant, std::forward<const std::vector<HypernodeID>>(partners));
    }
    if (_context.enable_soft_rebalance) {
      _soft_rebalance_refiner.postUncontraction(representant, std::forward<const std::vector<HypernodeID>>(partners));
    }
    if (_context.only_do_advanced_moves) {
      _local_search_refiner_am.postUncontraction(representant, std::forward<const std::vector<HypernodeID>>(partners));
    } else {
      _local_search_refiner.postUncontraction(representant, std::forward<const std::vector<HypernodeID>>(partners));
    }
  }

  void printSummary() const override {
    if (_context.partition.quiet_mode) {
      return;
    }

    if (_context.enable_soft_rebalance) {
      LOG << "Soft Rebalance Refiner:";
      _soft_rebalance_refiner.printSummary();
    } else {
      LOG << "Soft Rebalance disabled";
    }
    if (_context.enable_hard_rebalance) {
      LOG << "Hard Rebalance Refiner:";
      _hard_rebalance_refiner.printSummary();
    } else {
      LOG << "Hard Rebalance disabled";
    }
    if (_context.only_do_advanced_moves) {
      LOG << "AM Local Search Refiner:";
      _local_search_refiner_am.printSummary();
    } else {
      LOG << "GM Local Search Refiner:";
      _local_search_refiner.printSummary();
    }
    LOG << "Running with extra refinement nodes:" << _context.refine_rebalance_moves;
  }

private:
  void initializeImpl(const HyperedgeWeight max_gain) final {
    _qg->rebuild();
    _gain_manager->initialize();
    if (_context.only_do_advanced_moves) {
      _local_search_refiner_am.initialize(max_gain);
    } else {
      _local_search_refiner.initialize(max_gain);
    }
    if (_context.enable_hard_rebalance) {
      _hard_rebalance_refiner.initialize(max_gain);
    }
    if (_context.enable_soft_rebalance) {
      _soft_rebalance_refiner.initialize(max_gain);
    }
    _least_num_nodes = _hg.currentNumNodes();
    _is_initialized = true;

    // relax balance constrain until initial partition is initially within the bound
    _context.partition.epsilon = std::max(metrics::imbalance(_hg, _context), _context.partition.epsilon);
    _context.setupPartWeights(_hg.totalWeight());

    _context.only_do_advanced_moves = false;
    _context.repeated_insert = false;

    if (!_context.partition.quiet_mode) {
      if (_context.only_do_advanced_moves) {
        LOG << "\t-move heuristic: advanced moves";
      } else {
        if (_context.repeated_insert) {
          LOG << "\t-move heuristic: global moves with repeated insert";
        } else {
          LOG << "\t-move heuristic: global moves";
        }
      }
      LOG << "\t-extra refinement nodes:" << _context.refine_rebalance_moves;
      LOG << "\t-epsilon:" << _context.partition.epsilon;
    }
  }

  void performMovesAndUpdateCacheImpl(const std::vector<Move> &moves,
                                      std::vector<HypernodeID> &refinement_nodes,
                                      const UncontractionGainChanges &changes) final {
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

  bool refineImpl(std::vector<HypernodeID> &refinement_nodes,
                  const std::array<HypernodeWeight, 2> &max_allowed_part_weights,
                  const UncontractionGainChanges &uncontraction_changes,
                  Metrics &best_metrics) final {
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
      ASSERT(!_gain_manager->hasDelta());

      const auto &moves = _soft_rebalance_refiner.moves();
      if (_context.enable_hard_rebalance) {
        _hard_rebalance_refiner.performMovesAndUpdateCache(moves, refinement_nodes, uncontraction_changes);
      }
      if (_context.only_do_advanced_moves) {
        _local_search_refiner_am.performMovesAndUpdateCache(moves, refinement_nodes, uncontraction_changes);
      } else {
        _local_search_refiner.performMovesAndUpdateCache(moves, refinement_nodes, uncontraction_changes);
      }

      if (_context.refine_rebalance_moves) {
        for (const auto &move : moves) {
          if (std::find(more_refinement_nodes.begin(), more_refinement_nodes.end(), move.hn) ==
              more_refinement_nodes.end()) {
            more_refinement_nodes.push_back(move.hn);
          }
        }
      }
      DBG << "-->" << best_metrics.imbalance;
    }

    if (_context.enable_hard_rebalance && best_metrics.imbalance > _context.partition.epsilon) {
      DBG << "Running hard rebalance because" << best_metrics.imbalance << ">" << _context.partition.epsilon;
      _hard_rebalance_refiner.refine(refinement_nodes, max_allowed_part_weights, uncontraction_changes, best_metrics);
      ASSERT(!_gain_manager->hasDelta());

      const auto &moves = _hard_rebalance_refiner.moves();
      if (_context.enable_soft_rebalance) {
        _soft_rebalance_refiner.performMovesAndUpdateCache(moves, refinement_nodes, uncontraction_changes);
      }
      if (_context.only_do_advanced_moves) {
        _local_search_refiner_am.performMovesAndUpdateCache(moves, refinement_nodes, uncontraction_changes);
      } else {
        _local_search_refiner.performMovesAndUpdateCache(moves, refinement_nodes, uncontraction_changes);
      }

      if (_context.refine_rebalance_moves) {
        for (const auto &move : moves) {
          if (std::find(more_refinement_nodes.begin(), more_refinement_nodes.end(), move.hn) ==
              more_refinement_nodes.end()) {
            more_refinement_nodes.push_back(move.hn);
          }
        }
      }
      DBG << "-->" << best_metrics.imbalance;
    }

    const bool stop = _context.only_do_advanced_moves
                      ? _local_search_refiner_am.refine(more_refinement_nodes, max_allowed_part_weights,
                                                        uncontraction_changes, best_metrics)
                      : _local_search_refiner.refine(more_refinement_nodes, max_allowed_part_weights,
                                                     uncontraction_changes, best_metrics);
    ASSERT(!_gain_manager->hasDelta());

    const auto &moves = _context.only_do_advanced_moves
                        ? _local_search_refiner_am.moves()
                        : _local_search_refiner.moves();
    if (_context.enable_soft_rebalance) {
      _soft_rebalance_refiner.performMovesAndUpdateCache(moves, refinement_nodes, uncontraction_changes);
    }
    if (_context.enable_hard_rebalance) {
      _hard_rebalance_refiner.performMovesAndUpdateCache(moves, refinement_nodes, uncontraction_changes);
    }
    _qg->resetChanged();

    ASSERT(best_metrics.km1 == metrics::km1(_hg));
    ASSERT(best_metrics.imbalance == metrics::imbalance(_hg, _context));
    return stop;
  }

  Hypergraph &_hg;
  const Context &_original_context;
  Context _context;
  HypernodeID _least_num_nodes{0};
  std::shared_ptr<AdjacencyMatrixQuotientGraph<DFSCycleDetector>> _qg;
  std::shared_ptr<KMinusOneGainManager> _gain_manager;

  AcyclicLocalSearchRefiner<StoppingPolicy, FMImprovementPolicy> _local_search_refiner;
  //AcyclicLocalSearchRepeatedRefiner<StoppingPolicy, FMImprovementPolicy> _local_search_refiner;
  AcyclicKWayAdvancedMovesFMRefiner<StoppingPolicy, FMImprovementPolicy> _local_search_refiner_am;
  AcyclicHardRebalanceRefiner _hard_rebalance_refiner;
  AcyclicSoftRebalanceRefiner _soft_rebalance_refiner;
};

} // namespace kahypar