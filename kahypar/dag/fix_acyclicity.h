#pragma once

#include <stack>
#include <vector>
#include <utility>

#include "kahypar/dag/quotient_graph.h"
#include "kahypar/datastructure/hypergraph.h"
#include "kahypar/definitions.h"

namespace kahypar {
namespace dag {
namespace internal {
bool findCycleDFS(const CyclicQuotientGraph& cqg, QNodeID node, std::vector<std::pair<QNodeID, QNodeID>>& path, int mark, std::vector<int>& marked) {
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

void moveSuccessors(Hypergraph& hg, const HypernodeID hn, const PartitionID from_part, const PartitionID to_part) { // TODO can we improve on this?
  for (const HyperedgeID& he : hg.incidentTailEdges(hn)) {
    for (const HypernodeID& head : hg.heads(he)) {
      if (hg.partID(head) == to_part) {
        hg.changeNodePart(head, to_part, from_part);
        moveSuccessors(hg, head, from_part, to_part);
      }
    }
  }
}

void breakEdge(Hypergraph& hg, const PartitionID from_part, const PartitionID to_part) {
  std::vector<HypernodeID> hns_in_from_part;

  for (const HypernodeID& hn : hg.nodes()) {
    if (hg.partID(hn) == from_part) {
      hns_in_from_part.push_back(hn);
    }
  }

  for (const HypernodeID& hn : hns_in_from_part) {
    for (const HyperedgeID& he : hg.incidentTailEdges(hn)) {
      for (const HypernodeID& head : hg.heads(he)) {
        if (hg.partID(head) == to_part) {
          hg.changeNodePart(head, to_part, from_part);
          moveSuccessors(hg, head, from_part, to_part);
        }
      }
    }
  }

  ///////
  for (const HypernodeID& hn : hg.nodes()) {
    if (hg.partID(hn) != from_part) {
      continue;
    }

    for (const HyperedgeID& he : hg.incidentTailEdges(hn)) {
      for (const HypernodeID& head : hg.heads(he)) {
        if (hg.partID(head) == to_part) {
          LOG << hn << hg.partID(hn) << "-->" << head << hg.partID(head);
        }
      }
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
} // namespace internal

//void fixAcyclicity(Hypergraph& hg, const Context& context) {
//  CyclicQuotientGraph cqg(hg, context);
//
//  while (!cqg.isAcyclic()) {
//    const auto edge = internal::pickEdge(cqg);
//    internal::breakEdge(hg, cqg, edge.first, edge.second);
//  }
//
//  ASSERT(AdjacencyMatrixQuotientGraph<DFSCycleDetector>(hg, context).isAcyclic());
//}

void fixBipartitionAcyclicity(Hypergraph& hg, const Context& context, const PartitionID part0, const PartitionID part1) {
//  CyclicQuotientGraph cqg(hg, context);

//  while (!cqg.isAcyclic()) {
//    const auto edge = internal::pickEdge(cqg);
//    ASSERT(edge.first == part0 || edge.first == part1);
//    ASSERT(edge.second == part0 || edge.second == part1);
//    internal::breakEdge(hg, cqg, edge.first, edge.second);
//  }

  internal::breakEdge(hg, part0, part1);

  ASSERT(AdjacencyMatrixQuotientGraph<DFSCycleDetector>(hg, context).isAcyclic());
}

void fixBipartitionAcyclicity(Hypergraph &hg, const Context& context) {
  const auto original_partition = internal::createPartitionSnapshot(hg);

  // move successors of nodes in part 0 in part 1 to part 0
  internal::breakEdge(hg, 0, 1);
  double imbalance_01 = metrics::imbalance(hg, context);
  const auto partition_01 = internal::createPartitionSnapshot(hg);

  if (!AdjacencyMatrixQuotientGraph<DFSCycleDetector>(hg, context).isAcyclic()) {
    LOG << "Error, imbalance 01 is cyclic!";
    AdjacencyMatrixQuotientGraph<DFSCycleDetector>(hg, context).log();
    std::exit(1);
  }

  // other way around
  hg.setPartition(original_partition);
  internal::breakEdge(hg, 1, 0);
  double imbalance_10 = metrics::imbalance(hg, context);

  if (!AdjacencyMatrixQuotientGraph<DFSCycleDetector>(hg, context).isAcyclic()) {
    LOG << "Error, imbalance 10 is cyclic!";
    std::exit(1);
  }

  // select the most balanced partition from 01 and 10, i.e. if 01 is better balanced, restore it and otherwise do
  // nothing because we already have partition 10
  if (imbalance_01 <= imbalance_10) {
    LOG << "Chose 01";
    hg.setPartition(partition_01);
  } else {
    LOG << "Chose 10";
  }
}
} // namespace dag
} // namespace kabypar