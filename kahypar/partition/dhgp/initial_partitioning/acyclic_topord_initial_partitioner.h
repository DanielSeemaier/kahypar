#pragma once

#include <vector>

#include "kahypar/definitions.h"
#include "kahypar/partition/initial_partitioning/i_initial_partitioner.h"
#include "kahypar/partition/initial_partitioning/initial_partitioner_base.h"

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
  void partitionImpl() override final {
    Base::multipleRunsInitialPartitioning();
  }

  void initialPartition() {
  }

  using Base::_hg;
  using Base::_context;
};
} // namespace kahypar::dhgp