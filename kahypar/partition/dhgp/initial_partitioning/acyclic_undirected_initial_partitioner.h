#pragma once

#include <vector>

#include "kahypar/definitions.h"
#include "kahypar/partition/initial_partitioning/i_initial_partitioner.h"
#include "kahypar/partition/initial_partitioning/initial_partitioner_base.h"

namespace kahypar::dhgp {
class AcyclicUndirectedInitialPartitioner : public IInitialPartitioner,
                                 private InitialPartitionerBase<AcyclicUndirectedInitialPartitioner> {
  using Base = InitialPartitionerBase<AcyclicUndirectedInitialPartitioner>;
  friend Base;

 public:
  AcyclicUndirectedInitialPartitioner(Hypergraph& hypergraph, Context& context) :
      Base(hypergraph, context) { }

  ~AcyclicUndirectedInitialPartitioner() override = default;

  AcyclicUndirectedInitialPartitioner(const AcyclicUndirectedInitialPartitioner&) = delete;
  AcyclicUndirectedInitialPartitioner& operator= (const AcyclicUndirectedInitialPartitioner&) = delete;

  AcyclicUndirectedInitialPartitioner(AcyclicUndirectedInitialPartitioner&&) = delete;
  AcyclicUndirectedInitialPartitioner& operator= (AcyclicUndirectedInitialPartitioner&&) = delete;

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