/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2017 Sebastian Schlag <sebastian.schlag@kit.edu>
 * Copyright (C) 2017 Robin Andre <robinandre1995@web.de>
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

#include <algorithm>
#include <limits>
#include <utility>
#include <vector>

#include "kahypar/io/sql_plottools_serializer.h"
#include "kahypar/partition/evolutionary/edge_frequency.h"
#include "kahypar/partition/evolutionary/population.h"


namespace kahypar {
namespace combine {
static constexpr bool debug = false;

Individual partitions(Hypergraph& hg,
                      const Parents& parents,
                      Context& context) {
  const HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  DBG << V(context.evolutionary.action.decision());
  DBG << "Parent 1: initial" << V(parents.first.fitness());
  DBG << "Parent 2: initial" << V(parents.second.fitness());
  context.evolutionary.parent1 = &parents.first.partition();
  context.evolutionary.parent2 = &parents.second.partition();
#ifndef NDEBUG
  ASSERT(parents.first.fitness() == ([](Hypergraph& hg, const Parents& parents) -> int {
        hg.setPartition(parents.first.partition());
        HyperedgeWeight metric = metrics::km1(hg);
        hg.reset();
        return metric;
      })(hg,parents));
  DBG << "initial" << V(metrics::km1(hg)) << V(metrics::imbalance(hg, context));

  ASSERT(parents.second.fitness() == ([](Hypergraph& hg, const Parents& parents) -> int {
        hg.setPartition(parents.second.partition());
        HyperedgeWeight metric = metrics::km1(hg);
        hg.reset();
        return metric;
      })(hg,parents));

#endif

  hg.reset();
  const HypernodeID original_contraction_limit_multiplier =
    context.coarsening.contraction_limit_multiplier;
  if (context.evolutionary.unlimited_coarsening_contraction) {
    context.coarsening.contraction_limit_multiplier = 1;
  }

  Partitioner().partition(hg, context);

  const HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
  Timer::instance().add(context, Timepoint::evolutionary,
                        std::chrono::duration<double>(end - start).count());

  context.coarsening.contraction_limit_multiplier = original_contraction_limit_multiplier;
  DBG << "Offspring" << V(metrics::km1(hg)) << V(metrics::imbalance(hg, context));
  ASSERT(metrics::km1(hg) <= std::min(parents.first.fitness(), parents.second.fitness()));
  io::serializer::serializeEvolutionary(context, hg);
  return Individual(hg, context);
}


Individual usingTournamentSelection(Hypergraph& hg, const Context& context, const Population& population, const Population& imbalanced_population) {
  Context temporary_context(context);

  temporary_context.evolutionary.action =
    Action { meta::Int2Type<static_cast<int>(EvoDecision::combine)>() };
  temporary_context.coarsening.rating.rating_function = RatingFunction::heavy_edge;
  temporary_context.coarsening.rating.partition_policy = RatingPartitionPolicy::evolutionary;

  bool coin = Randomize::instance().flipCoin();
  if (coin && context.evolutionary.use_imbalanced_population) {
    if (imbalanced_population.size() == 0) {
      throw std::logic_error("imbalance_population is empty");
    }

    LOG << "Using tournament selection on population and imbalanced_population";
    const Parents parents = {
      population.singleTournamentSelection(),
      imbalanced_population.singleTournamentSelection()
    };
    return combine::partitions(hg, parents, temporary_context);
  } else if (coin && context.evolutionary.use_cross_combine) {
    LOG << "Using on-the-fly individual with random k and epsilon";
    const PartitionID lower_bound = std::max<PartitionID>(2, context.partition.k / 2);
    const PartitionID upper_bound = context.partition.k * 4;
    const PartitionID k_prime = Randomize::instance().getRandomInt(lower_bound, upper_bound);
    const double epsilon_prime = Randomize::instance().getRandomFloat(
      context.evolutionary.lower_epsilon_prime_bound,
      context.evolutionary.upper_epsilon_prime_bound
    );
    LOG << "\t" << V(k_prime) << V(epsilon_prime);

    Context individual_context(context);
    individual_context.partition.k = k_prime;
    individual_context.partition.epsilon = epsilon_prime;
    individual_context.partition.final_epsilon = epsilon_prime;
    individual_context.setupPartWeights(hg.totalWeight());

    hg.changeK(k_prime);
    Population individual_population;
    individual_population.generateIndividual(hg, individual_context);
    hg.changeK(context.partition.k);

    const Parents parents = {
      population.singleTournamentSelection(),
      individual_population.individualAt(0)
    };
    return combine::partitions(hg, parents, temporary_context);
  } else {
    LOG << "Using tournamentSelect() on population";
    const auto& parents = population.tournamentSelect();
    return combine::partitions(hg, parents, temporary_context);
  }
}



Individual edgeFrequency(Hypergraph& hg, const Context& context, const Population& population) {
  const HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  hg.reset();
  Context temporary_context(context);

  temporary_context.evolutionary.action =
    Action { meta::Int2Type<static_cast<int>(EvoDecision::combine)>(),
             meta::Int2Type<static_cast<int>(EvoCombineStrategy::edge_frequency)>() };

  temporary_context.coarsening.rating.rating_function = RatingFunction::edge_frequency;
  temporary_context.coarsening.rating.partition_policy = RatingPartitionPolicy::normal;
  temporary_context.coarsening.rating.heavy_node_penalty_policy =
    HeavyNodePenaltyPolicy::edge_frequency_penalty;

  temporary_context.evolutionary.edge_frequency =
    computeEdgeFrequency(population.listOfBest(context.evolutionary.edge_frequency_amount),
                         hg.initialNumEdges());

  DBG << V(temporary_context.evolutionary.action.decision());


  Partitioner().partition(hg, temporary_context);

  const HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
  Timer::instance().add(context, Timepoint::evolutionary,
                        std::chrono::duration<double>(end - start).count());


  DBG << "final result" << V(metrics::km1(hg)) << V(metrics::imbalance(hg, context));
  io::serializer::serializeEvolutionary(temporary_context, hg);
  return Individual(hg, context);
}
}  // namespace combine
}  // namespace kahypar
