#pragma once

#include <algorithm>

#include "kahypar/definitions.h"
#include "kahypar/partition/initial_partitioning/i_initial_partitioner.h"
#include "kahypar/partition/initial_partitioning/initial_partitioner_base.h"
#include "kahypar/utils/randomize.h"
#include "kahypar/dag/topological_ordering.h"

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
  void partitionImpl() final {
    Base::multipleRunsInitialPartitioning();
  }

  void initialPartition() {
    PartitionID unassigned_part = _context.initial_partitioning.unassigned_part;
    _context.initial_partitioning.unassigned_part = -1;
    Base::resetPartitioning();

    auto topological_ordering = dag::calculateTopologicalOrdering(_hg);

    // assign nodes to the lowest non-overloaded block in topological order
    PartitionID p = 0;
    for (const HypernodeID& hn : topological_ordering) {
      if (_hg.partWeight(p) >= _context.initial_partitioning.upper_allowed_partition_weight[p]) {
        ++p;
      }
      ASSERT(p < _context.partition.k, "Could not find a non-overloaded partition for node" << hn);
      _hg.setNodePart(hn, p);
    }
  }
};
} // namespace kahypar