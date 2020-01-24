/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2017 Sebastian Schlag <sebastian.schlag@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include <limits>
#include <vector>

#include "kahypar/definitions.h"
#include "kahypar/partition/coarsening/hypergraph_pruner.h"
#include "kahypar/partition/coarsening/i_coarsener.h"
#include "kahypar/partition/context.h"
#include "kahypar/partition/factories.h"
#include "kahypar/partition/metrics.h"
#include "kahypar/partition/multilevel.h"
#include "kahypar/partition/refinement/i_refiner.h"

namespace kahypar {
namespace acyclic_rb {
static constexpr bool debug = false;

static inline bool partitionVCycle(Hypergraph &hypergraph, ICoarsener &coarsener,
                                   IRefiner &refiner, const Context &context, const bool use_parent1 = false) {
  hypergraph.resetEdgeHashes();
  io::printVcycleBanner(context);
  io::printCoarseningBanner(context);

  coarsener.coarsen(context.coarsening.contraction_limit);
  if (context.partition.verbose_output && context.type == ContextType::main) {
    io::printHypergraphInfo(hypergraph, "Coarsened Hypergraph");
  }

  // special meaning of use_parent1 that doesnt make any sense, this is ugly, refactor this
  if ((context.partition_evolutionary && context.evolutionary.action.
    requires().evolutionary_parent_contraction) || use_parent1) {
    hypergraph.reset();
    hypergraph.setPartition(*context.evolutionary.parent1);

    if (!use_parent1) {
      const HyperedgeWeight parent_1_objective = metrics::correctMetric(hypergraph, context);

      hypergraph.setPartition(*context.evolutionary.parent2);
      const HyperedgeWeight parent_2_objective = metrics::correctMetric(hypergraph, context);

      if (parent_1_objective < parent_2_objective) {
        hypergraph.setPartition(*context.evolutionary.parent1);
      }
    }
  }

  hypergraph.initializeNumCutHyperedges();
  io::printLocalSearchBanner(context);
  const auto old_km1 = metrics::km1(hypergraph);
  const bool improved_quality = coarsener.uncoarsen(refiner);
  const auto new_km1 = metrics::km1(hypergraph);
  LOG << "vCycle result:" << old_km1 << "-->" << new_km1;
  io::printLocalSearchResults(context, hypergraph);
  return improved_quality;
}

static inline void performInitialPartitioning(Hypergraph &hypergraph, const Context &context) {
  recursive_bisection::partition(hypergraph, context);
}

static inline std::vector<PartitionID> createPartitionSnapshot(const Hypergraph &hg) {
  std::vector<PartitionID> partition(hg.initialNumNodes());
  for (const HypernodeID &hn : hg.nodes()) {
    partition[hn] = hg.partID(hn);
  }
  return partition;
}

static inline void partition(Hypergraph &hypergraph, const Context &context) {
  ASSERT(hypergraph.isDirected());

  Context ctx_copy(context);
  const bool recombine = context.partition_evolutionary && context.evolutionary.action.requires().initial_partitioning;
  std::vector<PartitionID> parent1;
  std::vector<PartitionID> parent2;

  if (recombine) {
    ctx_copy.evolutionary.action = Action{meta::Int2Type<static_cast<int>(EvoDecision::combine)>()};
    parent2 = createPartitionSnapshot(hypergraph);
    ctx_copy.evolutionary.parent2 = &parent2;
  }

  // perform initial partitioning
  if (!context.partition_evolutionary || context.evolutionary.action.requires().initial_partitioning) {
    if (context.initial_partitioning.level == InitialPartitioningLevel::finest) {
      hypergraph.resetPartitioning();
      performInitialPartitioning(hypergraph, context);
    } else {
      throw std::runtime_error("not implemented");
    }
  }

  if (recombine) {
    parent1 = createPartitionSnapshot(hypergraph);
    ctx_copy.evolutionary.parent1 = &parent1;
  }

  std::unique_ptr<ICoarsener> coarsener(
      CoarsenerFactory::getInstance().createObject(
          context.coarsening.algorithm, hypergraph, ctx_copy,
          hypergraph.weightOfHeaviestNode()));
  std::unique_ptr<IRefiner> km1_refiner(
      RefinerFactory::getInstance().createObject(
          context.local_search.algorithm, hypergraph, ctx_copy));

  if (metrics::imbalance(hypergraph, context) > context.partition.epsilon) {
    ctx_copy.partition.epsilon = metrics::imbalance(hypergraph, context);
    LOG << "Relaxed epsilon' from" << context.partition.epsilon << "to" << ctx_copy.partition.epsilon;
  }
  ctx_copy.setupPartWeights(hypergraph.totalWeight());

  LOG << "Initial" << context.partition.objective << "before first Vcycle ="
      << (context.partition.objective == Objective::cut
          ? metrics::hyperedgeCut(hypergraph)
          : metrics::km1(hypergraph));

  bool achieved_imbalance = metrics::imbalance(hypergraph, ctx_copy) < context.partition.final_epsilon;
  HyperedgeWeight best_km1 = metrics::km1(hypergraph);
  std::vector<PartitionID> best_partition = createPartitionSnapshot(hypergraph);

  for (uint32_t vcycle = 1; vcycle <= context.partition.global_search_iterations; ++vcycle) {
    LOG << "Performing vCycle on already partitioned graph using the following configuration:";
    LOG << "\t-refiner:" << context.local_search.algorithm;
    LOG << "\t-coarsening:" << context.coarsening.algorithm;
    LOG << "\t-recombine operation:" << recombine;
    LOG << "\t-advanced moves:" << context.only_do_advanced_moves;

    context.partition.current_v_cycle = vcycle;
    partitionVCycle(hypergraph, *coarsener, *km1_refiner, ctx_copy, recombine);

    const HyperedgeWeight current_km1 = metrics::km1(hypergraph);
    const bool current_achieved_imbalance = metrics::imbalance(hypergraph, ctx_copy) < context.partition.final_epsilon;

    // if previous partition is imbalanced, always accept the current partition if it is balanced or better than the previous one
    // otherwise only accept balanced partition that is better than the best balanced partition found so far
    if (best_km1 == 0 || (current_km1 < best_km1 && achieved_imbalance) ||
        (!achieved_imbalance && (current_achieved_imbalance || current_km1 < best_km1))) {
      LOG << "V-Cycle" << vcycle << "improved KM1 from" << best_km1 << "to" << current_km1;
      best_km1 = current_km1;
      best_partition = createPartitionSnapshot(hypergraph);
      achieved_imbalance = current_achieved_imbalance;
    } else {
      LOG << "V-Cycle" << vcycle << "did not improve KM1:" << best_km1 << "is better than" << current_km1
          << "or imbalance constrain is violated";
    }
  }

  // rollback to best partition found so far
  hypergraph.setPartition(best_partition);
  hypergraph.initializeNumCutHyperedges();

  // rebalance edge case: graph too small for multilevel
  const double imbalance = metrics::imbalance(hypergraph, context);
  if (imbalance > context.partition.final_epsilon || imbalance > context.partition.epsilon) {
    LOG << "Running hard rebalance to improve imbalance from" << imbalance << "to min{"
        << context.partition.final_epsilon << "," << context.partition.epsilon << "}";

    Context balanced_context = context;
    balanced_context.partition.epsilon = context.partition.final_epsilon;
    balanced_context.setupPartWeights(hypergraph.totalWeight());

    AdjacencyMatrixQuotientGraph<DFSCycleDetector> qg(hypergraph, balanced_context);
    KMinusOneGainManager gain_manager(hypergraph, balanced_context);
    gain_manager.initialize();
    AcyclicHardRebalanceRefiner hard_balance_refiner(hypergraph, balanced_context, qg, gain_manager);
    hard_balance_refiner.initialize(0);
    Metrics current_metrics = {metrics::hyperedgeCut(hypergraph),
                               metrics::km1(hypergraph),
                               metrics::imbalance(hypergraph, context)};
    UncontractionGainChanges changes{};
    std::vector<HypernodeID> refinement_nodes{};
    hard_balance_refiner.refine(refinement_nodes, {0, 0}, changes, current_metrics);

    hard_balance_refiner.printSummary();
  }
}
}  // namespace direct_kway
}  // namespace kahypar
