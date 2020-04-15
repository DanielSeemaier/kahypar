#pragma once

#ifndef KAHYPAR_ENABLE_DHGP
#error This file should not be included when building without KAHYPAR_ENABLE_DHGP
#endif // KAHYPAR_ENABLE_DHGP

#include <vector>

#include "kahypar/definitions.h"
#include "kahypar/partition/initial_partitioning/i_initial_partitioner.h"
#include "kahypar/partition/initial_partitioning/initial_partitioner_base.h"
#include "kahypar/partition/dhgp/topord.h"

namespace kahypar::dhgp {
class AcyclicTopOrdInitialPartitioner : public IInitialPartitioner,
                                 private InitialPartitionerBase<AcyclicTopOrdInitialPartitioner> {
  using Base = InitialPartitionerBase<AcyclicTopOrdInitialPartitioner>;
  friend Base;

 public:
  AcyclicTopOrdInitialPartitioner(Hypergraph& hypergraph, Context& context) :
    Base(hypergraph, context) { }

  ~AcyclicTopOrdInitialPartitioner() override = default;

  AcyclicTopOrdInitialPartitioner(const AcyclicTopOrdInitialPartitioner&) = delete;
  AcyclicTopOrdInitialPartitioner& operator= (const AcyclicTopOrdInitialPartitioner&) = delete;

  AcyclicTopOrdInitialPartitioner(AcyclicTopOrdInitialPartitioner&&) = delete;
  AcyclicTopOrdInitialPartitioner& operator= (AcyclicTopOrdInitialPartitioner&&) = delete;

 private:
  void partitionImpl() final {
    Base::multipleRunsInitialPartitioning();
  }

  void initialPartition() {
    const PartitionID unassigned_part = _context.initial_partitioning.unassigned_part;
    _context.initial_partitioning.unassigned_part = -1;
    Base::resetPartitioning();
    _context.initial_partitioning.unassigned_part = unassigned_part;

    const auto topord = calculateTopologicalOrdering(_hg, true);
    PartitionID p = 0;
    for (const HypernodeID& hn : topord) {
      if (_hg.partWeight(p) > _context.initial_partitioning.perfect_balance_partition_weight[p]) {
        ++p;
      }
      ASSERT(p < _context.partition.k);
      _hg.setNodePart(hn, p);
    }

    _hg.initializeNumCutHyperedges();
  }

  using Base::_hg;
  using Base::_context;
};
} // namespace kahypar::dhgp