#pragma once

#ifndef KAHYPAR_ENABLE_DHGP
#error This file should not be included when building without KAHYPAR_ENABLE_DHGP
#endif // KAHYPAR_ENABLE_DHGP

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
class AcyclicTwoWayFMRefiner final : public IRefiner,
                                     private FMRefinerBase<HypernodeID,
                                                           AcyclicTwoWayFMRefiner<StoppingPolicy,
                                                                                  FMImprovementPolicy>> {
 private:
  static constexpr bool enable_heavy_assert = false;
  static constexpr bool debug = false;

  using Base = FMRefinerBase<HypernodeID, AcyclicTwoWayFMRefiner<StoppingPolicy,
                                                                 FMImprovementPolicy>>;

  friend class FMRefinerBase<HypernodeID, AcyclicTwoWayFMRefiner<StoppingPolicy,
                                                                 FMImprovementPolicy>>;

  using HypernodeWeightArray = std::array<HypernodeWeight, 2>;
  using HEState = typename Base::HEState;
  using Base::kInvalidGain;
  using Base::kInvalidHN;

 public:
  AcyclicTwoWayFMRefiner(Hypergraph& hypergraph, const Context& context) : Base(hypergraph, context) {}

  ~AcyclicTwoWayFMRefiner() override = default;

  AcyclicTwoWayFMRefiner(const AcyclicTwoWayFMRefiner&) = delete;
  AcyclicTwoWayFMRefiner& operator=(const AcyclicTwoWayFMRefiner&) = delete;

  AcyclicTwoWayFMRefiner(AcyclicTwoWayFMRefiner&&) = delete;
  AcyclicTwoWayFMRefiner& operator=(AcyclicTwoWayFMRefiner&&) = delete;

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
      return false;
  }

  using Base::_hg;
  using Base::_context;
  using Base::_pq;
  using Base::_performed_moves;
  using Base::_hns_to_activate;
};
} // namespace kahypar::dhgp