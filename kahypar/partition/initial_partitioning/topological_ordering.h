#pragma once

#include <algorithm>
#include <kahypar/partition/refinement/acyclic_2way_km1_refiner.h>

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

    _hg.initializeNumCutHyperedges();

    if (_context.partition.k == 2) {
      Metrics current_metrics = {metrics::hyperedgeCut(_hg),
                                 metrics::km1(_hg),
                                 metrics::imbalance(_hg, _context)};
      HyperedgeWeight previous_km1 = current_metrics.km1;
      UncontractionGainChanges changes{}; // dummy

      do {
        previous_km1 = current_metrics.km1;

        AcyclicTwoWayKMinusOneRefiner local_search_refiner(_hg, _context);
        local_search_refiner.initialize(0);

        std::vector<HypernodeID> local_search_refinement_nodes;
        for (const HypernodeID& hn : _hg.nodes()) {
          if (_hg.isBorderNode(hn)) {
            local_search_refinement_nodes.push_back(hn);
          }
        }

        local_search_refiner.refine(local_search_refinement_nodes, {0, 0}, changes, current_metrics);
        //local_search_refiner.printSummary();
        LOG << "Result of 2Way refinement:" << previous_km1 << "-->" << current_metrics.km1;
      } while (0.99 * previous_km1 > current_metrics.km1);
    }
  }
};
} // namespace kahypar