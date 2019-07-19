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

#include "kahypar/dag/topological_ordering.h"
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
class MLCoarsener final : public ICoarsener,
                          private VertexPairCoarsenerBase<> {
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
  MLCoarsener(Hypergraph& hypergraph, const Context& context,
              const HypernodeWeight weight_of_heaviest_node) :
    Base(hypergraph, context, weight_of_heaviest_node),
    _rater(_hg, _context) {}

  ~MLCoarsener() override = default;

  MLCoarsener(const MLCoarsener&) = delete;
  MLCoarsener& operator=(const MLCoarsener&) = delete;

  MLCoarsener(MLCoarsener&&) = delete;
  MLCoarsener& operator=(MLCoarsener&&) = delete;

 private:
  void coarsenImpl(const HypernodeID limit) override final {
    int pass_nr = 0;
    const auto top = dag::calculateWeakTopologicalOrdering(_hg);

    std::vector<HypernodeID> current_hns;
    while (_hg.currentNumNodes() > limit) {
      DBG << V(pass_nr);
      DBG << V(_hg.currentNumNodes());
      DBG << V(_hg.currentNumEdges());
      _rater.resetMatches();
      current_hns.clear();
      const HypernodeID num_hns_before_pass = _hg.currentNumNodes();
      for (const HypernodeID& hn : _hg.nodes()) {
        current_hns.push_back(hn);
      }
      Randomize::instance().shuffleVector(current_hns, current_hns.size());

      // std::sort(std::begin(current_hns), std::end(current_hns),
      //           [&](const HypernodeID l, const HypernodeID r) {
      //             return _hg.nodeDegree(l) < _hg.nodeDegree(r);
      //           });

      // TODO move these modifications this to own coarsening algorithm
      std::vector<bool> marked(_hg.initialNumNodes(), false);

      for (const HypernodeID& u : current_hns) {
        if (marked[u] || !_hg.nodeIsEnabled(u)) {
          continue;
        }

        const Rating rating = _rater.rate(u);
        const HypernodeID v = rating.target;
        if (v == kInvalidTarget || marked[v]) {
          continue;
        }

        bool accept = top[u] + 1 == top[v]; // v --> u
        if (!accept) {
          bool v_is_succ = false;
          bool other_node_is_succ = false;
          for (const HyperedgeID& he : _hg.incidentTailEdges(u)) {
            for (const HypernodeID& head : _hg.heads(he)) {
              if (v == head) {
                v_is_succ = true;
              } else {
                other_node_is_succ = true;
              }
            }
          }
          accept = v_is_succ && !other_node_is_succ;
        }
        if (!accept) {
          bool u_is_pred = false;
          bool other_node_is_pred = false;
          for (const HyperedgeID& he : _hg.incidentHeadEdges(v)) {
            for (const HypernodeID& tail : _hg.tails(he)) {
              if (u == tail) {
                u_is_pred = true;
              } else {
                other_node_is_pred = true;
              }
            }
          }
          accept = u_is_pred && !other_node_is_pred;
        }

        if (!accept) {
          continue;
        }

        bool v_is_pred = false;
        for (const HyperedgeID& he : _hg.incidentHeadEdges(u)) {
          for (const HypernodeID& tail : _hg.tails(he)) {
            if (tail == v) {
              v_is_pred = true;
              break;
            }
          }
          if (v_is_pred) {
            break;
          }
        }

        _rater.markAsMatched(u);
        _rater.markAsMatched(rating.target);
        LOG << "Contracting" << V(u) << V(v);

        if (v_is_pred) {
          performContraction(v, u);
          if (!dag::isAcyclic(_hg)) { LOG << "Error:" << V(u) << V(v); std::exit(1); }

          for (const HyperedgeID& he : _hg.incidentTailEdges(v)) {
            for (const HypernodeID& head : _hg.heads(he)) {
              if (top[v] == top[head] - 1) {
                marked[head] = false;
              }
            }
          }
        } else {
          performContraction(u, v);
          if (!dag::isAcyclic(_hg)) { LOG << "Error:" << V(u) << V(v); std::exit(1); }

          for (const HyperedgeID& he : _hg.incidentTailEdges(u)) {
            for (const HypernodeID& head : _hg.heads(he)) {
              if (top[u] == top[head] - 1) {
                marked[head] = false;
              }
            }
          }
        }

        marked[u] = true;
        marked[v] = true;

        if (_hg.currentNumNodes() <= limit) {
          break;
        }
      }

      if (num_hns_before_pass == _hg.currentNumNodes()) {
        break;
      }
      ++pass_nr;
    }
    _context.stats.add(StatTag::Coarsening, "HnsAfterCoarsening", _hg.currentNumNodes());
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
