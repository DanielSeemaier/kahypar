/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2016 Sebastian Schlag <sebastian.schlag@kit.edu>
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
#include <string>
#include <vector>

#include "kahypar/definitions.h"
#include "kahypar/macros.h"
#include "kahypar/partition/coarsening/policies/fixed_vertex_acceptance_policy.h"
#include "kahypar/partition/coarsening/policies/rating_acceptance_policy.h"
#include "kahypar/partition/coarsening/policies/rating_community_policy.h"
#include "kahypar/partition/coarsening/policies/rating_heavy_node_penalty_policy.h"
#include "kahypar/partition/coarsening/policies/rating_partition_policy.h"
#include "kahypar/partition/coarsening/policies/rating_score_policy.h"
#include "kahypar/partition/coarsening/policies/rating_tie_breaking_policy.h"
#include "kahypar/partition/coarsening/vertex_pair_rater.h"
#include "kahypar/dag/acyclic_clustering.h"
#include "kahypar/utils/htimer.h"

namespace kahypar {
template <class ScorePolicy = HeavyEdgeScore,
    class HeavyNodePenaltyPolicy = NoWeightPenalty,
    class CommunityPolicy = UseCommunityStructure,
    class RatingPartitionPolicy = NormalPartitionPolicy,
    class AcceptancePolicy = BestRatingPreferringUnmatched<>,
    class FixedVertexPolicy = AllowFreeOnFixedFreeOnFreeFixedOnFixed,
    typename RatingType = RatingType>
class AcyclicClusterCoarseningV2 final : public ICoarsener, private VertexPairCoarsenerBase<>{
private:
  static constexpr bool debug = false;

  static constexpr HypernodeID kInvalidTarget = std::numeric_limits<HypernodeID>::max();

  using Base = VertexPairCoarsenerBase;
  using Rater = VertexPairRater<ScorePolicy,
      HeavyNodePenaltyPolicy,
      CommunityPolicy,
      RatingPartitionPolicy,
      AcceptancePolicy,
      FixedVertexPolicy,
      RatingType>;
  using Rating = typename Rater::Rating;

public:
  AcyclicClusterCoarseningV2(Hypergraph& hypergraph, const Context& context,
                             const HypernodeWeight weight_of_heaviest_node) :
      Base(hypergraph, context, weight_of_heaviest_node),
      _rater(_hg, _context) { }

  ~AcyclicClusterCoarseningV2() override = default;

  AcyclicClusterCoarseningV2(const AcyclicClusterCoarseningV2&) = delete;
  AcyclicClusterCoarseningV2& operator= (const AcyclicClusterCoarseningV2&) = delete;

  AcyclicClusterCoarseningV2(AcyclicClusterCoarseningV2&&) = delete;
  AcyclicClusterCoarseningV2& operator= (AcyclicClusterCoarseningV2&&) = delete;

private:
  void coarsenImpl(const HypernodeID limit) override final {
    std::size_t num_contractions = 0;
    HTimer timer;
    timer.start();
    while (_hg.currentNumNodes() > limit) {
      num_contractions = 0;

      const double max_weight_fraction = static_cast<double>(_context.coarsening.max_allowed_weight_multiplier)
          / static_cast<double>(_context.coarsening.contraction_limit_multiplier * _context.partition.k);
      const auto clustering = dag::findAcyclicClusteringWithCycleDetection(_hg, _context, max_weight_fraction);
      for (const HypernodeID& hn : _hg.nodes()) {
        if (clustering[hn] != hn) {
          performContraction(clustering[hn], hn);
          ++num_contractions;
        }
      }

//      if (!dag::isAcyclic(_hg)) {
//        throw std::runtime_error("coarser graph is cyclic");
//      } else {
//        LOG << "Iteration 1 OK, got" << _hg.currentNumNodes() << "/" << limit << "/" << _hg.initialNumNodes();
//      }

      if (num_contractions == 0) {
        break;
      }
    }

    if (!dag::isAcyclic(_hg)) {
      throw std::runtime_error("coarser graph is cyclic");
    }

    _context.stats.add(StatTag::Coarsening, "HnsAfterCoarsening", _hg.currentNumNodes());
    LOG << "Coarsened from" << _hg.initialNumNodes() << "to" << _hg.currentNumNodes() << "in" << timer.stop();
  }

  bool uncoarsenImpl(IRefiner& refiner) override final {
    return doUncoarsen(refiner);
  }

  using Base::_pq;
  using Base::_hg;
  using Base::_context;
  using Base::_history;
  Rater _rater;
};
}  // namespace kahypar