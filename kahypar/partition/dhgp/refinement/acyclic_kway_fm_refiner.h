#pragma once

#ifndef KAHYPAR_ENABLE_DHGP
#error This file should not be included when building without KAHYPAR_ENABLE_DHGP
#endif// KAHYPAR_ENABLE_DHGP

#include <algorithm>
#include <array>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "kahypar/datastructure/binary_heap.h"
#include "kahypar/datastructure/fast_reset_array.h"
#include "kahypar/datastructure/fast_reset_flag_array.h"
#include "kahypar/datastructure/sparse_set.h"
#include "kahypar/definitions.h"
#include "kahypar/meta/mandatory.h"
#include "kahypar/meta/template_parameter_to_string.h"
#include "kahypar/partition/context.h"
#include "kahypar/partition/metrics.h"
#include "kahypar/partition/refinement/2way_fm_gain_cache.h"
#include "kahypar/partition/refinement/fm_refiner_base.h"
#include "kahypar/partition/refinement/i_refiner.h"
#include "kahypar/partition/refinement/move.h"
#include "kahypar/partition/refinement/policies/fm_improvement_policy.h"
#include "kahypar/utils/float_compare.h"
#include "kahypar/utils/randomize.h"

namespace kahypar::dhgp {
template<class StoppingPolicy = Mandatory,
         class FMImprovementPolicy = CutDecreasedOrInfeasibleImbalanceDecreased>
class AcyclicKWayFMRefiner final : public IRefiner,
                                   private FMRefinerBase<HypernodeID,
                                                         AcyclicKWayFMRefiner<StoppingPolicy,
                                                                              FMImprovementPolicy>> {
 private:
  static constexpr bool enable_heavy_assert = false;
  static constexpr bool debug = false;

  using Base = FMRefinerBase<HypernodeID, AcyclicKWayFMRefiner<StoppingPolicy,
                                                               FMImprovementPolicy>>;

  friend class FMRefinerBase<HypernodeID, AcyclicKWayFMRefiner<StoppingPolicy,
                                                               FMImprovementPolicy>>;

  using HEState = typename Base::HEState;
  using Base::kInvalidGain;
  using Base::kInvalidHN;

 public:
  AcyclicKWayFMRefiner(Hypergraph& hypergraph, const Context& context) : Base(hypergraph, context) {}

  ~AcyclicKWayFMRefiner() override = default;

  AcyclicKWayFMRefiner(const AcyclicKWayFMRefiner&) = delete;
  AcyclicKWayFMRefiner& operator=(const AcyclicKWayFMRefiner&) = delete;

  AcyclicKWayFMRefiner(AcyclicKWayFMRefiner&&) = delete;
  AcyclicKWayFMRefiner& operator=(AcyclicKWayFMRefiner&&) = delete;

  bool isInitialized() const {
    return _is_initialized;
  }

 private:
  void initializeImpl(const HyperedgeWeight max_gain) override final {
  }

  void performMovesAndUpdateCacheImpl(const std::vector<Move>& moves,
                                      std::vector<HypernodeID>& refinement_nodes,
                                      const UncontractionGainChanges& uncontraction_changes) override final {
  }

  std::vector<Move> rollbackImpl() override final {
    return std::vector<Move>();
  }

  bool refineImpl(std::vector<HypernodeID>& refinement_nodes,
                  const HypernodeWeightArray& max_allowed_part_weights,
                  const UncontractionGainChanges& changes,
                  Metrics& best_metrics) override final {
  }

  using Base::_context;
  using Base::_hg;
  using Base::_hns_to_activate;
  using Base::_performed_moves;
  using Base::_pq;
};
}// namespace kahypar::dhgp