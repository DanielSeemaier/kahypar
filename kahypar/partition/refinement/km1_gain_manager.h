#pragma once

#include "kahypar/datastructure/sparse_map.h"
#include "kahypar/datastructure/fast_reset_array.h"
#include "kway_fm_gain_cache.h"
#include "kahypar/datastructure/hypergraph.h"
#include "kahypar/partition/context.h"

namespace kahypar {
class KMinusOneGainManager {
 private:
  static constexpr bool debug = false;
  static constexpr HypernodeID hn_to_debug = 2;

  using GainCache = KwayGainCache<Gain>;

  struct PinState {
    char one_pin_in_from_part_before : 1;
    char one_pin_in_to_part_after : 1;
    char two_pins_in_from_part_before : 1;
    char two_pins_in_to_part_after : 1;

    PinState(const bool one_in_from_before, const bool one_in_to_after,
             const bool two_in_from_before, const bool two_in_to_after) :
      one_pin_in_from_part_before(one_in_from_before),
      one_pin_in_to_part_after(one_in_to_after),
      two_pins_in_from_part_before(two_in_from_before),
      two_pins_in_to_part_after(two_in_to_after) {}
  };

  struct NonadjacentGainCacheDelta {
    const HypernodeID hn;
    const Gain delta;

    NonadjacentGainCacheDelta(const HypernodeID hn, const Gain delta) :
      hn(hn),
      delta(delta) {}
  };

 public:
  KMinusOneGainManager(const Hypergraph& hg, const Context& context) :
    _hg(hg),
    _context(context),
    _adjacent_gain_cache(_hg.initialNumNodes(), _context.partition.k),
    _nonadjacent_gain_cache(_hg.initialNumNodes()),
    _tmp_gains(_context.partition.k, 0),
    _new_adjacent_part(_hg.initialNumNodes(), Hypergraph::kInvalidPartition) {}

  KMinusOneGainManager(const KMinusOneGainManager&) = delete;
  KMinusOneGainManager& operator=(const KMinusOneGainManager&) = delete;

  KMinusOneGainManager(KMinusOneGainManager&&) = delete;
  KMinusOneGainManager& operator=(KMinusOneGainManager&&) = delete;

  ~KMinusOneGainManager() = default;

  void initialize() {
    _adjacent_gain_cache.clear();
    for (const HypernodeID& hn : _hg.nodes()) {
      initializeGainCacheFor(hn);
    }
  }

  void preUncontraction(const HypernodeID representant) {
    _adjacent_gain_cache.clear(representant);
    resetDelta();
  }

  void postUncontraction(const HypernodeID representant, const std::vector<HypernodeID>& partners) {
    initializeGainCacheFor(representant);
    for (const HypernodeID partner : partners) {
      ASSERT(!_adjacent_gain_cache.entryExists(partner), V(representant) << V(partner));
      initializeGainCacheFor(partner);
    }
    resetDelta();
  }

  void updateAfterMovement(const HypernodeID moved_hn, const PartitionID from_part, const PartitionID to_part) {
    updateAfterMovement(moved_hn, from_part, to_part,
                        [](const auto&, const auto&, const auto&) {},
                        [](const auto&, const auto&, const auto&) {},
                        [](const auto&, const auto&, const auto&) {});
  }

  template<typename ConnectivityIncreaseCallback, typename UpdateCallback, typename ConnectivityDecreaseCallback>
  void updateAfterMovement(const HypernodeID moved_hn, const PartitionID from_part, const PartitionID to_part,
                           ConnectivityIncreaseCallback&& connectivity_inc_cb, UpdateCallback&& update_cb,
                           ConnectivityDecreaseCallback&& connectivity_dec_cb) {
    DBGC(moved_hn == hn_to_debug) << "updateAfterMovement(" << moved_hn << "," << from_part << "," << to_part << ")";
    ASSERT(_hg.partID(moved_hn) == to_part);
    ASSERT(!_adjacent_gain_cache.entryExists(moved_hn, from_part));
    _new_adjacent_part.resetUsedEntries();
    const Gain old_nonadjacent_gain = _nonadjacent_gain_cache[moved_hn];

    const bool moved_hn_internal_before = isInternal(moved_hn);
    bool moved_hn_remains_connected_to_from_part = false;

    for (const HyperedgeID& he : _hg.incidentEdges(moved_hn)) {
      const HypernodeID pins_in_source_part_after = _hg.pinCountInPart(he, from_part);
      const HypernodeID pins_in_target_part_after = _hg.pinCountInPart(he, to_part);
      const HyperedgeWeight he_weight = _hg.edgeWeight(he);

      moved_hn_remains_connected_to_from_part |= pins_in_source_part_after != 0;

      if (pins_in_source_part_after != 0 && pins_in_target_part_after == 1) {
        _nonadjacent_gain_cache[moved_hn] += he_weight;
        _nonadjacent_gain_cache_deltas.emplace_back(moved_hn, he_weight);
        update_cb(moved_hn, Hypergraph::kInvalidPartition, he_weight);
      } else if (pins_in_source_part_after == 0 && pins_in_target_part_after != 1) {
        _nonadjacent_gain_cache[moved_hn] -= he_weight;
        _nonadjacent_gain_cache_deltas.emplace_back(moved_hn, -he_weight);
        update_cb(moved_hn, Hypergraph::kInvalidPartition, -he_weight);
      }

      if (!moved_hn_internal_before) {
        if (pins_in_source_part_after == 0 && _hg.pinCountInPart(he, to_part) != 1) {
          for (const PartitionID& part : _adjacent_gain_cache.adjacentParts(moved_hn)) {
            if (part != from_part && part != to_part) {
              _adjacent_gain_cache.updateExistingEntry(moved_hn, part, -_hg.edgeWeight(he));
              update_cb(moved_hn, part, -_hg.edgeWeight(he));
            }
          }
        } else if (pins_in_source_part_after != 0 && _hg.pinCountInPart(he, to_part) == 1) {
          for (const PartitionID& part : _adjacent_gain_cache.adjacentParts(moved_hn)) {
            if (part != from_part && part != to_part) {
              _adjacent_gain_cache.updateExistingEntry(moved_hn, part, _hg.edgeWeight(he));
              update_cb(moved_hn, part, _hg.edgeWeight(he));
            }
          }
        }
      }

      for (const HypernodeID& pin : _hg.pins(he)) {
        if (pin == moved_hn) {
          continue;
        }

        const bool move_decreased_connectivity = (_hg.pinCountInPart(he, from_part) == 0);
        const bool move_increased_connectivity = (_hg.pinCountInPart(he, to_part) == 1);
        performConnectivityUpdate(pin, from_part, to_part, move_decreased_connectivity, move_increased_connectivity,
                                  connectivity_inc_cb, connectivity_dec_cb);
        performDeltaGainUpdate(pin, from_part, to_part, he, update_cb);
      }
    }

    if (moved_hn_internal_before && moved_hn_remains_connected_to_from_part) {
      DBGC(moved_hn == hn_to_debug)
      << "Moved HN" << moved_hn << "was previously internal, but after moving it from" << from_part << "to" << to_part
      << "remains connected to" << from_part << "and is thus inserted into the adjacent gain cache";


      _adjacent_gain_cache.addEntryDueToConnectivityIncrease(moved_hn, from_part, -old_nonadjacent_gain);
    } else {
      if (!_adjacent_gain_cache.entryExists(moved_hn, to_part)) {
        _adjacent_gain_cache.addEntryDueToConnectivityIncrease(moved_hn, to_part, old_nonadjacent_gain);
      }

      ASSERT(_adjacent_gain_cache.entryExists(moved_hn, to_part));
      _adjacent_gain_cache.updateFromAndToPartOfMovedHN(moved_hn, from_part, to_part,
                                                        moved_hn_remains_connected_to_from_part);
    }
  }

  Gain gain(const HypernodeID hn, const PartitionID to_part) const {
    ASSERT(_hg.partID(hn) != to_part);

    if (_adjacent_gain_cache.entryExists(hn, to_part)) {
      return adjacentGain(hn, to_part);
    }

    return nonadjacentGain(hn);
  }

  Gain adjacentGain(const HypernodeID hn, const PartitionID to_part) const {
    ASSERT(_hg.partID(hn) != to_part);
    return _adjacent_gain_cache.entry(hn, to_part);
  }

  Gain nonadjacentGain(const HypernodeID hn) const {
    return _nonadjacent_gain_cache[hn];
  }

  const auto& adjacentParts(const HypernodeID hn) const {
    ASSERT(_adjacent_gain_cache.entryExists(hn));
    return _adjacent_gain_cache.adjacentParts(hn);
  }

  bool isInternal(const HypernodeID hn) const {
    return !_adjacent_gain_cache.nonemptyEntryExists(hn);
  }

  bool isAdjacentTo(const HypernodeID hn, const PartitionID part) const {
    ASSERT(_hg.partID(hn) != part);
    const bool adjacent_to = _adjacent_gain_cache.entryExists(hn, part);
    ASSERT(adjacent_to == isConnectedTo(hn, part));
    return adjacent_to;
  }

  void resetDelta() {
    _adjacent_gain_cache.resetDelta();
    _nonadjacent_gain_cache_deltas.clear();
  }

  void rollbackDelta() {
    _adjacent_gain_cache.rollbackDelta();
    for (auto rit = _nonadjacent_gain_cache_deltas.crbegin(); rit != _nonadjacent_gain_cache_deltas.crend(); ++rit) {
      _nonadjacent_gain_cache[rit->hn] -= rit->delta;
    }
    _nonadjacent_gain_cache_deltas.clear();
  }

 private:
  template<typename ConnectivityIncreaseCallback, typename ConnectivityDecreaseCallback>
  void performConnectivityUpdate(const HypernodeID pin, const PartitionID from_part, const PartitionID to_part,
                                 const bool move_decreased_connectivity, const bool move_increased_connectivity,
                                 ConnectivityIncreaseCallback&& connectivity_inc_cb,
                                 ConnectivityDecreaseCallback&& connectivity_dec_cb) {
    if (move_decreased_connectivity && _adjacent_gain_cache.entryExists(pin, from_part) &&
        !isConnectedTo(pin, from_part)) {
      _adjacent_gain_cache.removeEntryDueToConnectivityDecrease(pin, from_part);
      connectivity_dec_cb(pin, from_part, nonadjacentGain(pin));
    }
    if (move_increased_connectivity && !_adjacent_gain_cache.entryExists(pin, to_part)) {
      const Gain gain = gainInducedByHypergraph(pin, to_part);
      _adjacent_gain_cache.addEntryDueToConnectivityIncrease(pin, to_part, gain);
      _new_adjacent_part.set(pin, to_part);
      connectivity_inc_cb(pin, to_part, gain);
    }
  }

  template<typename UpdateCallback>
  void performDeltaGainUpdate(const HypernodeID pin, const PartitionID from_part, const PartitionID to_part,
                              const HyperedgeID he, UpdateCallback&& update_cb) {
    const HypernodeID pin_count_from_part_before_move = _hg.pinCountInPart(he, from_part) + 1;
    const HypernodeID pin_count_to_part_after_move = _hg.pinCountInPart(he, to_part);
    const bool move_decreased_connectivity = pin_count_from_part_before_move - 1 == 0;
    const bool move_increased_connectivity = pin_count_to_part_after_move == 1;
    const HyperedgeWeight he_weight = _hg.edgeWeight(he);

    const PinState pin_state(pin_count_from_part_before_move == 1,
                             pin_count_to_part_after_move == 1,
                             pin_count_from_part_before_move == 2,
                             pin_count_to_part_after_move == 2);

    const PartitionID source_part = _hg.partID(pin);

    if (source_part == from_part) {
      if (pin_state.two_pins_in_from_part_before) {
        _nonadjacent_gain_cache[pin] += he_weight;
        _nonadjacent_gain_cache_deltas.emplace_back(pin, he_weight);
        update_cb(pin, Hypergraph::kInvalidPartition, he_weight);

        for (const PartitionID& part : _adjacent_gain_cache.adjacentParts(pin)) {
          if (_new_adjacent_part.get(pin) != part) {
            _adjacent_gain_cache.updateExistingEntry(pin, part, he_weight);
            update_cb(pin, part, he_weight);
          }
        }
      }
    } else if (source_part == to_part && pin_state.two_pins_in_to_part_after) {
      _nonadjacent_gain_cache[pin] -= he_weight;
      _nonadjacent_gain_cache_deltas.emplace_back(pin, -he_weight);
      update_cb(pin, Hypergraph::kInvalidPartition, -he_weight);

      for (const PartitionID& part : _adjacent_gain_cache.adjacentParts(pin)) {
        if (_new_adjacent_part.get(pin) != part) {
          _adjacent_gain_cache.updateExistingEntry(pin, part, -he_weight);
          update_cb(pin, part, -he_weight);
        }
      }
    }

    if (pin_state.one_pin_in_from_part_before && _adjacent_gain_cache.entryExists(pin, from_part)) {
      _adjacent_gain_cache.updateExistingEntry(pin, from_part, -he_weight);
      update_cb(pin, from_part, -he_weight);
    }

    if (pin_state.one_pin_in_to_part_after && _new_adjacent_part.get(pin) != to_part) {
      _adjacent_gain_cache.updateExistingEntry(pin, to_part, he_weight);
      update_cb(pin, to_part, he_weight);
    }
  }

  void initializeGainCacheFor(const HypernodeID hn) {
    _tmp_gains.clear();
    const PartitionID source_part = _hg.partID(hn);
    HyperedgeWeight internal = 0;
    for (const HyperedgeID& he : _hg.incidentEdges(hn)) {
      const HyperedgeWeight he_weight = _hg.edgeWeight(he);
      internal += _hg.pinCountInPart(he, source_part) != 1 ? he_weight : 0;
      for (const PartitionID& part : _hg.connectivitySet(he)) {
        _tmp_gains[part] += he_weight;
      }
    }

    for (const auto& target_part : _tmp_gains) {
      if (target_part.key == source_part) {
        ASSERT(!_adjacent_gain_cache.entryExists(hn, source_part));
        continue;
      }

      _adjacent_gain_cache.initializeEntry(hn, target_part.key, target_part.value - internal);
    }

    _nonadjacent_gain_cache_deltas.emplace_back(hn, -_nonadjacent_gain_cache[hn] - internal);
    _nonadjacent_gain_cache[hn] = -internal;
  }

  bool isConnectedTo(const HypernodeID hn, const PartitionID part) const {
    for (const HyperedgeID& he : _hg.incidentEdges(hn)) {
      if (_hg.pinCountInPart(he, part) > 0) {
        return true;
      }
    }

    return false;
  }

  Gain gainInducedByHypergraph(const HypernodeID hn, const PartitionID target_part) const {
    Gain gain = 0;
    for (const HyperedgeID& he : _hg.incidentEdges(hn)) {
      ASSERT(_hg.edgeSize(he) > 1, V(he));
      gain += gainInducedByHyperedge(hn, he, target_part);
    }
    return gain;
  }

  Gain gainInducedByHyperedge(const HypernodeID hn, const HyperedgeID he, const PartitionID target_part) const {
    const HypernodeID pins_in_source_part = _hg.pinCountInPart(he, _hg.partID(hn));
    const HypernodeID pins_in_target_part = _hg.pinCountInPart(he, target_part);
    const HyperedgeWeight he_weight = _hg.edgeWeight(he);
    Gain gain = pins_in_source_part == 1 ? he_weight : 0;
    gain -= pins_in_target_part == 0 ? he_weight : 0;
    return gain;
  }

  const Hypergraph& _hg;
  const Context& _context;

  GainCache _adjacent_gain_cache;
  std::vector<Gain> _nonadjacent_gain_cache;
  ds::SparseMap<PartitionID, Gain> _tmp_gains;
  ds::FastResetArray<PartitionID> _new_adjacent_part;
  std::vector<NonadjacentGainCacheDelta> _nonadjacent_gain_cache_deltas;
};
}

