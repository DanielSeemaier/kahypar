#pragma once

#include <stack>
#include <vector>
#include <utility>

#include "kahypar/partition/refinement/acyclic_2way_km1_refiner.h"
#include "kahypar/partition/refinement/wip_refiner.h"
#include "kahypar/dag/quotient_graph.h"
#include "kahypar/datastructure/hypergraph.h"
#include "kahypar/definitions.h"

namespace kahypar {
namespace dag {
namespace internal {
bool
findCycleDFS(const CyclicQuotientGraph& cqg, QNodeID node, std::vector<std::pair<QNodeID, QNodeID>>& path, int mark,
             std::vector<int>& marked) {
  for (const QNodeID& child : cqg.outs(node)) {
    path.emplace_back(node, child);
    if (marked[child] == mark) { // detected cycle
      return true;
    } else {
      marked[child] = mark;
      if (findCycleDFS(cqg, child, path, mark, marked)) {
        return true;
      }
      path.pop_back();
    }
  }

  return false;
}

std::vector<std::pair<QNodeID, QNodeID>> findCycle(const CyclicQuotientGraph& cqg) {
  std::vector<std::pair<QNodeID, QNodeID>> cycle;
  std::vector<std::pair<QNodeID, QNodeID>> path;
  std::vector<int> marked(cqg.numberOfNodes());
  int mark = 1;
  for (QNodeID node = 0; node < cqg.numberOfNodes(); ++node) {
    if (marked[node] == 0) {
      if (findCycleDFS(cqg, node, path, mark, marked)) {
        break;
      }
      ++mark;
    }
  }
  ASSERT(!path.empty());

  const QNodeID witness = path.back().second;
  for (auto rit = path.crbegin(); rit != path.crend(); ++rit) {
    cycle.push_back(*rit);
    if (rit->first == witness) {
      break;
    }
  }
  std::reverse(cycle.begin(), cycle.end());

  ASSERT(cycle.size() > 1, V(cycle.size()));
  ASSERT([&]() {
    const auto& first = cycle.front();
    const auto& last = cycle.back();

    for (std::size_t i = 0; i < cycle.size(); ++i) {
      if (cqg.edgeWeight(cycle[i].first, cycle[i].second) <= 0) {
        return false;
      }

      if (i == 0) {
        if (cycle[i].first != last.second) {
          return false;
        }
      } else if (i + 1 == cycle.size()) {
        if (cycle[i].second != first.first) {
          return false;
        }
      } else {
        if (cycle[i].first != cycle[i - 1].second) {
          return false;
        }
      }
    }
    return true;
  }());

  return cycle;
}

std::pair<QNodeID, QNodeID> pickEdge(const CyclicQuotientGraph& cqg) {
  const auto cycle = findCycle(cqg);
  ASSERT(!cycle.empty());

  LLOG << cycle.front().first;
  for (const auto& edge : cycle) {
    LLOG << edge.second;
  }
  LOG << " ";

  std::pair<QNodeID, QNodeID> lightest_edge = cycle.front();
  for (const auto& edge : cycle) {
    if (cqg.edgeWeight(edge.first, edge.second) < cqg.edgeWeight(lightest_edge.first, lightest_edge.second)) {
      lightest_edge = edge;
    }
  }

  ASSERT([&]() {
    for (const auto& edge : cycle) {
      if (cqg.edgeWeight(edge.first, edge.second) < cqg.edgeWeight(lightest_edge.first, lightest_edge.second)) {
        return false;
      }
    }
    return true;
  }());
  return lightest_edge;
}

void moveSuccessors(Hypergraph& hg, const HypernodeID hn, const PartitionID from_part,
                    const PartitionID to_part) {
  for (const HyperedgeID& he : hg.incidentTailEdges(hn)) {
    for (const HypernodeID& head : hg.heads(he)) {
      if (hg.partID(head) == from_part) {
        hg.changeNodePart(head, from_part, to_part);
        moveSuccessors(hg, head, from_part, to_part);
      }
    }
  }
}

void movePredecessors(Hypergraph &hg, const HypernodeID hn, const PartitionID from_part,
        const PartitionID to_part) {
    for (const HyperedgeID& he : hg.incidentHeadEdges(hn)) {
        for (const HypernodeID& tail : hg.tails(he)) {
            if (hg.partID(tail) == from_part) {
                hg.changeNodePart(tail, from_part, to_part);
                movePredecessors(hg, tail, from_part, to_part);
            }
        }
    }
}

void breakEdge(Hypergraph& hg, const PartitionID u, const PartitionID v, const bool direction) {
  if (direction) { // move successors from v to u
      for (const HypernodeID &hn : hg.nodes()) {
          if (hg.partID(hn) != u) {
              continue;
          }
          moveSuccessors(hg, hn, v, u); // move all successors from v to u
      }
  } else { // move predecessors from u to v
      for (const HypernodeID& hn : hg.nodes()) {
          if (hg.partID(hn) != v) {
              continue;
          }
          movePredecessors(hg, hn, u, v);
      }
  }
}

std::vector<PartitionID> createPartitionSnapshot(const Hypergraph& hg) {
  std::vector<PartitionID> partition;
  partition.reserve(hg.currentNumNodes());
  for (const HypernodeID& hn : hg.nodes()) {
    partition.push_back(hg.partID(hn));
  }
  return partition;
}

void rebalancePartition(Hypergraph& hg, const Context& context, const bool refine_km1 = false) {
  hg.initializeNumCutHyperedges();

  Context balanced_context = context;
  balanced_context.partition.epsilon = context.partition.final_epsilon;
  balanced_context.setupPartWeights(hg.totalWeight());

  Metrics current_metrics = {metrics::hyperedgeCut(hg),
                             metrics::km1(hg),
                             metrics::imbalance(hg, context)};
  AdjacencyMatrixQuotientGraph<DFSCycleDetector> qg(hg, balanced_context);
  KMinusOneGainManager gain_manager(hg, balanced_context);
  gain_manager.initialize();
  UncontractionGainChanges changes{}; // dummy

  if (current_metrics.imbalance > balanced_context.partition.epsilon) {
    AcyclicHardRebalanceRefiner hard_balance_refiner(hg, balanced_context, qg, gain_manager);
    hard_balance_refiner.initialize(0);
    std::vector<HypernodeID> refinement_nodes{};
    hard_balance_refiner.refine(refinement_nodes, {0, 0}, changes, current_metrics);
  }

  if (refine_km1) {
    HyperedgeWeight previous_km1 = current_metrics.km1;
    do {
      previous_km1 = current_metrics.km1;

//      AcyclicLocalSearchRepeatedRefiner<NumberOfFruitlessMovesStopsSearch> local_search_refiner(hg, context);
      AcyclicTwoWayKMinusOneRefiner local_search_refiner(hg, context);
      local_search_refiner.initialize(0);

      std::vector<HypernodeID> local_search_refinement_nodes;
      for (const HypernodeID& hn : hg.nodes()) {
        if (hg.isBorderNode(hn)) {
          local_search_refinement_nodes.push_back(hn);
        }
      }

      local_search_refiner.refine(local_search_refinement_nodes, {0, 0}, changes, current_metrics);
      if (!context.partition.quiet_mode) {
        LOG << "Result of flat 2Way refinement:" << previous_km1 << "-->" << current_metrics.km1;
      }
    } while (0.99 * previous_km1 > current_metrics.km1);
  }
}
} // namespace internal

void fixBipartitionAcyclicity(Hypergraph& hg, const Context& context) {
  const auto original_partition = internal::createPartitionSnapshot(hg);
  const auto km_initial = metrics::km1(hg);

  // variant 1
  internal::breakEdge(hg, 0, 1, false);
  const auto imbalance_010_prime = metrics::imbalance(hg, context);
  const auto km_010_prime = metrics::km1(hg);
  if (context.initial_partitioning.balance_partition) {
    internal::rebalancePartition(hg, context, true);
  }
  double imbalance_010 = metrics::imbalance(hg, context);
  const auto km_010 = metrics::km1(hg);
  const auto partition_010 = internal::createPartitionSnapshot(hg);
  ASSERT(AdjacencyMatrixQuotientGraph<DFSCycleDetector>(hg, context).isAcyclic());

  // variant 2
  hg.setPartition(original_partition);
  internal::breakEdge(hg, 0, 1, true);
  const auto imbalance_011_prime = metrics::imbalance(hg, context);
  const auto km_011_prime = metrics::km1(hg);
  if (context.initial_partitioning.balance_partition) {
    internal::rebalancePartition(hg, context, true);
  }
  double imbalance_011 = metrics::imbalance(hg, context);
  const auto km_011 = metrics::km1(hg);
  const auto partition_011 = internal::createPartitionSnapshot(hg);
  ASSERT(AdjacencyMatrixQuotientGraph<DFSCycleDetector>(hg, context).isAcyclic());

  // variant 3
  hg.setPartition(original_partition);
  internal::breakEdge(hg, 1, 0, false);
  const auto imbalance_100_prime = metrics::imbalance(hg, context);
  const auto km_100_prime = metrics::km1(hg);
  if (context.initial_partitioning.balance_partition) {
    internal::rebalancePartition(hg, context, true);
  }
  double imbalance_100 = metrics::imbalance(hg, context);
  const auto km_100 = metrics::km1(hg);
  const auto partition_100 = internal::createPartitionSnapshot(hg);
  ASSERT(AdjacencyMatrixQuotientGraph<DFSCycleDetector>(hg, context).isAcyclic());

  // variant 4
  hg.setPartition(original_partition);
  internal::breakEdge(hg, 1, 0, true);
  const auto imbalance_101_prime = metrics::imbalance(hg, context);
  const auto km_101_prime = metrics::km1(hg);
  if (context.initial_partitioning.balance_partition) {
    internal::rebalancePartition(hg, context, true);
  }
  double imbalance_101 = metrics::imbalance(hg, context);
  const auto km_101 = metrics::km1(hg);
  ASSERT(AdjacencyMatrixQuotientGraph<DFSCycleDetector>(hg, context).isAcyclic());

//  LOG << V(km_initial) << V(km_010) << V(imbalance_010) << V(km_010_prime) << V(imbalance_010_prime);
//  LOG << V(km_initial) << V(km_011) << V(imbalance_011) << V(km_011_prime) << V(imbalance_011_prime);
//  LOG << V(km_initial) << V(km_100) << V(imbalance_100) << V(km_100_prime) << V(imbalance_100_prime);
//  LOG << V(km_initial) << V(km_101) << V(imbalance_101) << V(km_101_prime) << V(imbalance_101_prime);

  // apply best partition
  const auto best_km1 = std::min(km_010, std::min(km_011, std::min(km_100, km_101)));
  if (best_km1 == km_101) {
      // already set
  } else if (best_km1 == km_100) {
      hg.setPartition(partition_100);
  } else if (best_km1 == km_011) {
      hg.setPartition(partition_011);
  } else if (best_km1 == km_010) {
      hg.setPartition(partition_010);
  }
}
} // namespace dag
} // namespace kabypar