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

#include <fstream>
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

namespace kahypar {
template<class ScorePolicy = HeavyEdgeScore,
    class HeavyNodePenaltyPolicy = NoWeightPenalty,
    class CommunityPolicy = UseCommunityStructure,
    class RatingPartitionPolicy = NormalPartitionPolicy,
    class AcceptancePolicy = BestRatingPreferringUnmatched<>,
    class FixedVertexPolicy = AllowFreeOnFixedFreeOnFreeFixedOnFixed,
    typename RatingType = RatingType>
class ExternalCoarsener final : public ICoarsener,
                                private VertexPairCoarsenerBase<> {
private:
  static constexpr bool debug = true;

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
  ExternalCoarsener(Hypergraph &hypergraph, const Context &context,
                    const HypernodeWeight weight_of_heaviest_node) :
      Base(hypergraph, context, weight_of_heaviest_node),
      _rater(_hg, _context) {}

  ~ExternalCoarsener() override = default;

  ExternalCoarsener(const ExternalCoarsener &) = delete;

  ExternalCoarsener &operator=(const ExternalCoarsener &) = delete;

  ExternalCoarsener(ExternalCoarsener &&) = delete;

  ExternalCoarsener &operator=(ExternalCoarsener &&) = delete;

private:
  void coarsenImpl(const HypernodeID limit) override final {
    std::ifstream file(_context.coarsening.external_file);
    if (!file) {
      throw std::runtime_error("External coarsening file does not exist");
    }

    int u = 0, v = 0;
    LOG << "Using matching from" << _context.coarsening.external_file;
    LOG << "\t" << 0 << _hg.initialNumNodes();
    while (file >> u) {
      if (u < 0) { // separator
        const HypernodeID iteration = _hg.initialNumNodes() - _hg.currentNumNodes();
        _coarsening_levels.push_back(iteration);
        LOG << "\t" << _coarsening_levels.size() << _hg.currentNumNodes();
        continue;
      }

      file >> v;

      if (!_hg.nodeIsEnabled(u)) {
        throw std::runtime_error("Rep node is not enabled");
      }
      if (!_hg.nodeIsEnabled(v)) {
        throw std::runtime_error("Contract node is not enabled");
      }

      performContraction(
          static_cast<HypernodeID>(u),
          static_cast<HypernodeID>(v)
      );
    }

    LOG << "Contracted" << _hg.initialNumNodes() - _hg.currentNumNodes() << "nodes";
    _context.stats.add(StatTag::Coarsening, "HnsAfterCoarsening", _hg.currentNumNodes());
  }

  bool uncoarsenImpl(IRefiner &refiner) override final {
    return doUncoarsen(refiner);
  }

  using Base::_pq;
  using Base::_hg;
  using Base::_context;
  using Base::_history;
  Rater _rater;
  std::vector<HypernodeID> _coarsening_levels;
};
}  // namespace kahypar