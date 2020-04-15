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

namespace kahypar::dhgp {
template<class ScorePolicy = HeavyEdgeScore,
         class HeavyNodePenaltyPolicy = NoWeightPenalty,
         class CommunityPolicy = UseCommunityStructure,
         class RatingPartitionPolicy = NormalPartitionPolicy,
         class AcceptancePolicy = BestRatingPreferringUnmatched<>,
         class FixedVertexPolicy = AllowFreeOnFixedFreeOnFreeFixedOnFixed,
         typename RatingType = RatingType>
class AcyclicCoarsener final : public ICoarsener,
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
  AcyclicCoarsener(Hypergraph& hypergraph, const Context& context,
                   const HypernodeWeight weight_of_heaviest_node) : Base(hypergraph, context, weight_of_heaviest_node),
                                                                    _rater(_hg, _context) {}

  ~AcyclicCoarsener() override = default;

  AcyclicCoarsener(const AcyclicCoarsener&) = delete;
  AcyclicCoarsener& operator=(const AcyclicCoarsener&) = delete;

  AcyclicCoarsener(AcyclicCoarsener&&) = delete;
  AcyclicCoarsener& operator=(AcyclicCoarsener&&) = delete;

 private:
  void coarsenImpl(const HypernodeID limit) override final {
  }

  bool uncoarsenImpl(IRefiner& refiner) override final {
    return doUncoarsen(refiner);
  }

  using Base::_context;
  using Base::_hg;
  using Base::_history;
  using Base::_pq;
  Rater _rater;
};
}// namespace kahypar::dhgp
