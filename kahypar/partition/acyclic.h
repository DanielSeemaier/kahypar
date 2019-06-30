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
  return improved_quality;
}

static inline void performInitialPartitioning(Hypergraph& hypergraph, const Context& context, IRefiner& refiner) {
  io::printInitialPartitioningBanner(context);

  auto start = std::chrono::high_resolution_clock::now();
  initial::partition(hypergraph, context);
  LOG << "Initial" << context.partition.objective << " before initial refinement ="
      << (context.partition.objective == Objective::cut
          ? metrics::hyperedgeCut(hypergraph)
          : metrics::km1(hypergraph));
  LOG << "Performing initial refinement";

  hypergraph.initializeNumCutHyperedges();
  refiner.initialize(0);
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

  refiner.refine(refinement_nodes, {0, 0}, changes, current_metrics);
  ASSERT(QuotientGraph<DFSCycleDetector>(hypergraph, context).isAcyclic(), "Initial partition is not acyclic!");

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

static inline void partition(Hypergraph& hypergraph, const Context& context) {
  std::unique_ptr<ICoarsener> coarsener(
    CoarsenerFactory::getInstance().createObject(
      context.coarsening.algorithm, hypergraph, context,
      hypergraph.weightOfHeaviestNode()));
  std::unique_ptr<IRefiner> km1_refiner(
    RefinerFactory::getInstance().createObject(
      context.local_search.algorithm, hypergraph, context));
  std::unique_ptr<IRefiner> ip_refiner(
    RefinerFactory::getInstance().createObject(
      context.initial_partitioning.local_search.algorithm, hypergraph, context));

  if (!context.partition.vcycle_refinement_for_input_partition) {
    performInitialPartitioning(hypergraph, context, *ip_refiner);
  }

  for (uint32_t vcycle = 1; vcycle <= context.partition.global_search_iterations; ++vcycle) {
    context.partition.current_v_cycle = vcycle;
    const bool improved_quality = partitionVCycle(hypergraph, *coarsener, *km1_refiner, context);
    ASSERT(QuotientGraph<DFSCycleDetector>(hypergraph, context).isAcyclic(),
           "Vcycle" << vcycle << "produced a cyclic partition");

    if (!improved_quality) {
      LOG << "No improvement in V-cycle" << vcycle << ". Stopping global search.";
      break;
    }
  }

  // rebalance edge case: graph too small for multilevel
  const double imbalance = metrics::imbalance(hypergraph, context);
  if (imbalance > context.partition.final_epsilon || imbalance > context.partition.epsilon) {
    LOG << "Running hard rebalance to improve imbalance from" << imbalance << "to min{"
        << context.partition.final_epsilon << "," << context.partition.epsilon << "}";

    AcyclicHardRebalanceRefiner hard_balance_refiner(hypergraph, context);
    UncontractionGainChanges changes;
    changes.representative.push_back(0);
    changes.contraction_partner.push_back(0);
    hard_balance_refiner.initialize(0);
    Metrics current_metrics = {metrics::hyperedgeCut(hypergraph),
                               metrics::km1(hypergraph),
                               metrics::imbalance(hypergraph, context)};
    std::vector<HypernodeID> refinement_nodes{};
    hard_balance_refiner.changeEpsilon(context.partition.final_epsilon);
    hard_balance_refiner.refine(refinement_nodes, {0, 0}, changes, current_metrics);
  }
}
}  // namespace direct_kway
}  // namespace kahypar
