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
namespace acyclic {
static constexpr bool debug = false;

static inline bool partitionVCycle(Hypergraph& hypergraph, ICoarsener& coarsener,
                                   IRefiner& refiner, const Context& context) {
  // In order to perform parallel net detection, we have to reset the edge hashes
  // before coarsening.
  hypergraph.resetEdgeHashes();

  io::printVcycleBanner(context);
  io::printCoarseningBanner(context);

  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  coarsener.coarsen(context.coarsening.contraction_limit);
  HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
  Timer::instance().add(context, Timepoint::v_cycle_coarsening,
                        std::chrono::duration<double>(end - start).count());

  if (context.partition.verbose_output && context.type == ContextType::main) {
    io::printHypergraphInfo(hypergraph, "Coarsened Hypergraph");
  }

  hypergraph.initializeNumCutHyperedges();

  io::printLocalSearchBanner(context);

  start = std::chrono::high_resolution_clock::now();
  const bool improved_quality = coarsener.uncoarsen(refiner);
  end = std::chrono::high_resolution_clock::now();
  Timer::instance().add(context, Timepoint::v_cycle_local_search,
                        std::chrono::duration<double>(end - start).count());

  io::printLocalSearchResults(context, hypergraph);
  refiner.printSummary();
  return improved_quality;
}

static inline void performInitialPartitioning(Hypergraph& hypergraph, const Context& context) {
  io::printInitialPartitioningBanner(context);

  auto start = std::chrono::high_resolution_clock::now();
  initial::partition(hypergraph, context);
  LOG << "Initial" << context.partition.objective << " before initial refinement ="
      << (context.partition.objective == Objective::cut
          ? metrics::hyperedgeCut(hypergraph)
          : metrics::km1(hypergraph));

  hypergraph.initializeNumCutHyperedges();

  Context ctx_copy = context;
  std::unique_ptr<IRefiner> ip_refiner(
    RefinerFactory::getInstance().createObject(
      context.initial_partitioning.local_search.algorithm, hypergraph, ctx_copy));
  if (metrics::imbalance(hypergraph, context) > context.partition.epsilon) {
    ctx_copy.partition.epsilon = metrics::imbalance(hypergraph, context);
  }
  ctx_copy.setupPartWeights(hypergraph.totalWeight());

  if (context.partition.refine_initial_partition) {
    LOG << "Performing initial refinement with configured refiner";
    ip_refiner->initialize(0);
    UncontractionGainChanges changes;
    changes.representative.push_back(0);
    changes.contraction_partner.push_back(0);
    std::vector<HypernodeID> refinement_nodes;
    for (const HypernodeID& hn : hypergraph.nodes()) {
      refinement_nodes.push_back(hn);
    }
    Metrics current_metrics = {metrics::hyperedgeCut(hypergraph),
                               metrics::km1(hypergraph),
                               metrics::imbalance(hypergraph, context)};
    ip_refiner->refine(refinement_nodes, {0, 0}, changes, current_metrics);
    ASSERT(AdjacencyMatrixQuotientGraph<DFSCycleDetector>(hypergraph, context).isAcyclic(),
           "Initial partition is not acyclic!");
  } else {
    LOG << "Initial refinement disabled";
  }

  auto end = std::chrono::high_resolution_clock::now();
  Timer::instance().add(context, Timepoint::initial_partitioning,
                        std::chrono::duration<double>(end - start).count());

  if (context.partition.verbose_output && context.type == ContextType::main) {
    LOG << "Initial Partitioning Result:";
    LOG << "Initial" << context.partition.objective << "      ="
        << (context.partition.objective == Objective::cut
            ? metrics::hyperedgeCut(hypergraph)
            : metrics::km1(hypergraph));
    LOG << "Initial imbalance =" << metrics::imbalance(hypergraph, context);
    LOG << "Initial part sizes and weights:";
    io::printPartSizesAndWeights(hypergraph);
    LLOG << "Target weights:";
    LLOG << "w(*) =" << context.partition.max_part_weights[0] << "\n";
  }
}

static inline std::vector<PartitionID> createPartitionSnapshot(const Hypergraph& hg) {
  std::vector<PartitionID> partition(hg.initialNumNodes());
  for (const HypernodeID& hn : hg.nodes()) {
    partition[hn] = hg.partID(hn);
  }
  return partition;
}

static inline void partition(Hypergraph& hypergraph, const Context& context) {
  ASSERT(hypergraph.isDirected());

  Context ctx_copy = context;
  ctx_copy.enable_hard_rebalance = context.enable_hard_rebalance;

  std::unique_ptr<ICoarsener> coarsener(
    CoarsenerFactory::getInstance().createObject(
      context.coarsening.algorithm, hypergraph, context,
      hypergraph.weightOfHeaviestNode()));
  std::unique_ptr<IRefiner> km1_refiner(
    RefinerFactory::getInstance().createObject(
      context.local_search.algorithm, hypergraph, ctx_copy));
  std::unique_ptr<IRefiner> km1_refiner_second(
    RefinerFactory::getInstance().createObject(
      context.local_search.algorithm_second, hypergraph, ctx_copy));

  if (!context.partition.vcycle_refinement_for_input_partition) {
    performInitialPartitioning(hypergraph, context);
  }

  if (metrics::imbalance(hypergraph, context) > context.partition.epsilon) {
    ctx_copy.partition.epsilon = metrics::imbalance(hypergraph, context);
    LOG << "Relaxed epsilon' from" << context.partition.epsilon << "to" << ctx_copy.partition.epsilon;
  }
  ctx_copy.setupPartWeights(hypergraph.totalWeight());

  LOG << "Initial" << context.partition.objective << " before first Vcycle ="
      << (context.partition.objective == Objective::cut
          ? metrics::hyperedgeCut(hypergraph)
          : metrics::km1(hypergraph));

  bool achieved_imbalance = metrics::imbalance(hypergraph, ctx_copy) < context.partition.final_epsilon;
  HyperedgeWeight best_km1 = metrics::km1(hypergraph);
  std::vector<PartitionID> best_partition = createPartitionSnapshot(hypergraph);

  for (uint32_t vcycle = 1; vcycle <= context.partition.global_search_iterations; ++vcycle) {
    context.partition.current_v_cycle = vcycle;
    if (vcycle == 1) {
      LOG << "Using Refiner:" << context.local_search.algorithm;
      partitionVCycle(hypergraph, *coarsener, *km1_refiner, ctx_copy);
      ASSERT(AdjacencyMatrixQuotientGraph<DFSCycleDetector>(hypergraph, context).isAcyclic(),
             "Vcycle" << vcycle << "produced a cyclic partition");
    } else {
      LOG << "Using Refiner:" << context.local_search.algorithm_second;
      ctx_copy.partition.epsilon = context.partition.epsilon;
      ctx_copy.setupPartWeights(hypergraph.totalWeight());
      LOG << "Reset epsilon' to" << context.partition.epsilon << "for" << vcycle << "-nd/th vcycle";

      partitionVCycle(hypergraph, *coarsener, *km1_refiner_second, context);
      ASSERT(AdjacencyMatrixQuotientGraph<DFSCycleDetector>(hypergraph, context).isAcyclic(),
             "Vcycle" << vcycle << "produced a cyclic partition");
    }

    const HyperedgeWeight current_km1 = metrics::km1(hypergraph);
    const bool current_achieved_imbalance = metrics::imbalance(hypergraph, ctx_copy) < context.partition.final_epsilon;
    if (current_km1 < best_km1 && (!achieved_imbalance || current_achieved_imbalance)) {
      LOG << "V-Cycle" << vcycle << "improved KM1 from" << best_km1 << "to" << current_km1;
      best_km1 = current_km1;
      best_partition = createPartitionSnapshot(hypergraph);
      achieved_imbalance = current_achieved_imbalance;
    } else {
      LOG << "V-Cycle" << vcycle << "did not improve KM1:" << best_km1 << "is better than" << current_km1 << "or imbalance constrain is violated";
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
