#pragma once

#include <algorithm>

#include "kahypar/definitions.h"
#include "kahypar/partition/initial_partitioning/i_initial_partitioner.h"
#include "kahypar/partition/initial_partitioning/initial_partitioner_base.h"
#include "kahypar/utils/randomize.h"

namespace kahypar {
class TopologicalOrderingInitialPartitioner : public IInitialPartitioner,
                                              private InitialPartitionerBase<TopologicalOrderingInitialPartitioner>{
  using Base = InitialPartitionerBase<TopologicalOrderingInitialPartitioner>;
  friend Base;

 public:
  TopologicalOrderingInitialPartitioner(Hypergraph& hypergraph, Context& context) :
    Base(hypergraph, context) { }

  ~TopologicalOrderingInitialPartitioner() override = default;

  TopologicalOrderingInitialPartitioner(const TopologicalOrderingInitialPartitioner&) = delete;
  TopologicalOrderingInitialPartitioner& operator= (const TopologicalOrderingInitialPartitioner&) = delete;

  TopologicalOrderingInitialPartitioner(TopologicalOrderingInitialPartitioner&&) = delete;
  TopologicalOrderingInitialPartitioner& operator= (TopologicalOrderingInitialPartitioner&&) = delete;

 private:
  void partitionImpl() override final {
    Base::multipleRunsInitialPartitioning();
  }

  void initialPartition() {
    PartitionID unassigned_part = _context.initial_partitioning.unassigned_part;
    _context.initial_partitioning.unassigned_part = -1;
    Base::resetPartitioning();

    std::vector<HypernodeID> position(_hg.currentNumNodes());
    for (const HyperedgeID &he : _hg.edges()) {
      for (const HypernodeID &hn : _hg.heads(he)) {
        position[hn] += _hg.edgeNumTails(he);
      }
    }

    std::vector<HypernodeID> candidates;
    for (const HypernodeID &hn : _hg.nodes()) {
      if (position[hn] == 0) {
        candidates.push_back(hn);
      }
    }

    HypernodeID target_block_size = _hg.currentNumNodes() / _context.partition.k;
    HypernodeID current_block_size = 0;
    PartitionID current_block = 0;

    while (!candidates.empty()) {
      // select random candidate and remove it from the candidates list
      HypernodeID u = Randomize::instance().popRandomElement(candidates);

      // assign node to the current block
      if (current_block_size > target_block_size) {
        current_block = std::min(_context.partition.k - 1, current_block + 1);
        current_block_size = 0;
      }
      _hg.setNodePart(u, current_block);
      ++current_block_size;

      // decrease position of heads of edges that contain u as tail
      for (const HyperedgeID &he : _hg.incidentEdges(u)) {
        if (!_hg.isTail(u, he)) continue;

        for (const HypernodeID &hn : _hg.heads(he)) {
          ASSERT(position[hn] > 0);
          --position[hn];
          if (position[hn] == 0) {
            candidates.push_back(hn);
          }
        }
      }
    }

    ASSERT([&]() {
      for (const HypernodeID& hn : _hg.nodes()) {
        if (_hg.partID(hn) == -1) {
          return false;
        }
      }
      return true;
    } (), "There are unassigned hypernodes -- the graph contains cycles!");
  }
};
} // namespace kahypar