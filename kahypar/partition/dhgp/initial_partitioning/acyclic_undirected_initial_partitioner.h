#pragma once

#include <vector>

#include "kahypar/definitions.h"
#include "kahypar/utils/randomize.h"
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
    runUndirectedHgp();
    fixCyclicBipartition();
    _hg.initializeNumCutHyperedges();
  }

  void runUndirectedHgp() {
    const auto undirected_ctx = createContextForUndirectedPartitioner();
  }

  Context createContextForUndirectedPartitioner() const {
    Context undirected_ctx = _context;
    undirected_ctx.partition.mode = Mode::direct_kway;
    undirected_ctx.partition.seed = Randomize::instance().newRandomSeed();
    undirected_ctx.preprocessing.enable_min_hash_sparsifier = true;
    undirected_ctx.preprocessing.enable_community_detection = true;
    undirected_ctx.coarsening.algorithm = CoarseningAlgorithm::ml_style;
    undirected_ctx.initial_partitioning.mode = Mode::recursive_bisection;
    undirected_ctx.initial_partitioning.algo = InitialPartitionerAlgorithm::pool;
    undirected_ctx.initial_partitioning.local_search.algorithm = RefinementAlgorithm::twoway_fm;
    undirected_ctx.local_search.algorithm = RefinementAlgorithm::twoway_fm_hyperflow_cutter;
    return undirected_ctx;
  }

  void fixCyclicBipartition() {
    const auto cyclic_part = createPartitionSnapshot();
    using PartitionMetricPair = std::pair<std::vector<PartitionID>, HyperedgeWeight>;
    std::vector<PartitionMetricPair> acyclic_part;

    for (const PartitionID u : {0, 1}) {
      const PartitionID v = 1 - u;
      for (const bool direction : {false, true}) {
        breakQuotientGraphEdge(u, v, direction);
        // TODO balance, refine
        acyclic_part.emplace_back(createPartitionSnapshot(), metrics::km1(_hg));
        _hg.setPartition(cyclic_part);
      }
    }

    const auto best_part = std::min_element(acyclic_part.begin(), acyclic_part.end(), [](const auto& a, const auto& b) {
      return a.second < b.second;
    })->first;
    _hg.setPartition(best_part);
  }

  void breakQuotientGraphEdge(const PartitionID u, const PartitionID v, const bool direction) {
    if (direction) {
      for (const HypernodeID hn : _hg.nodes()) {
        if (_hg.partID(hn) == u) {
          moveSuccessors(hn, v, u);
        }
      }
    } else {
      for (const HypernodeID hn : _hg.nodes()) {
        if (_hg.partID(hn) == v) {
          movePredecessors(hn, u, v);
        }
      }
    }
  }

  void moveSuccessors(const HypernodeID hn, const PartitionID from, const PartitionID to) {
    for (const HyperedgeID& he : _hg.incidentTailEdges(hn)) {
      for (const HypernodeID& head : _hg.headPins(he)) {
        if (_hg.partID(head) == from) {
          _hg.changeNodePart(head, from, to);
          moveSuccessors(head, from, to);
        }
      }
    }
  }

  void movePredecessors(const HypernodeID hn, const PartitionID from, const PartitionID to) {
    for (const HyperedgeID& he : _hg.incidentHeadEdges(hn)) {
      for (const HypernodeID& tail : _hg.tailPins(he)) {
        if (_hg.partID(tail) == from) {
          _hg.changeNodePart(tail, from, to);
          movePredecessors(tail, from, to);
        }
      }
    }
  }

  std::vector<PartitionID> createPartitionSnapshot() const {
    std::vector<PartitionID> part;
    part.reserve(_hg.currentNumNodes());
    for (const HypernodeID hn : _hg.nodes()) {
      part.push_back(hn);
    }
    return part;
  }

  using Base::_hg;
  using Base::_context;
};
} // namespace kahypar::dhgp