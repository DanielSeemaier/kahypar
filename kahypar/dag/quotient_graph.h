#pragma once

#include <algorithm>
#include <functional>
#include <limits>
#include <memory>
#include <numeric>
#include <set>
#include <utility>
#include <vector>
#include <stack>

#include "gtest/gtest_prod.h"

#include "kahypar/macros.h"

#include "kahypar/datastructure/sparse_map.h"
#include "kahypar/definitions.h"
#include "kahypar/partition/context.h"
#include "kahypar/utils/randomize.h"
#include "kahypar/dag/topological_ordering.h"

namespace kahypar {
template<typename Detector>
class QuotientGraph {
 public:
  QuotientGraph(const Hypergraph& hypergraph, const Context& context) :
    _hg(hypergraph),
    _context(context),
    _adjacency_matrix(context.partition.k),
    _detector(context.partition.k) {
    constructFromHypergraph();
  }

  QuotientGraph(const QuotientGraph& other) = delete;
  QuotientGraph& operator=(const QuotientGraph& other) = delete;

  QuotientGraph(QuotientGraph&& other) noexcept = default;
  QuotientGraph& operator=(QuotientGraph&& other) = delete;

  ~QuotientGraph() = default;

  QNodeID numberOfNodes() const {
    return _adjacency_matrix.size();
  }

  QEdgeID numberOfEdges() const {
    return std::accumulate(_adjacency_matrix.begin(), _adjacency_matrix.end(), 0,
                           [](const QEdgeID& acc, const std::unordered_map<QNodeID, QEdgeWeight>& map) {
                             return acc + map.size();
                           });
  }

  bool update(const HypernodeID hn, const PartitionID from, const PartitionID to) {
    AdjacencyMatrix deltas(numberOfNodes());

    // find changes to edge weights of the quotient graph caused by this movement
    for (const HyperedgeID& he : _hg.incidentEdges(hn)) {
      if (_hg.isHead(hn, he)) { // hn is head
        for (const HypernodeID& ht : _hg.tails(he)) {
          if (_hg.partID(ht) != from) {
            deltas[_hg.partID(ht)][from] -=_hg.edgeWeight(he);
          }
          if (_hg.partID(ht) != to) {
            deltas[_hg.partID(ht)][to] += _hg.edgeWeight(he);
          }
        }
      } else { // hn is tail
        ASSERT(_hg.isTail(hn, he));
        for (const HypernodeID& hh : _hg.heads(he)) {
          if (_hg.partID(hh) != from) {
            deltas[from][_hg.partID(hh)] -= _hg.edgeWeight(he);
          }
          if (_hg.partID(hh) != to) {
            deltas[to][_hg.partID(hh)] += _hg.edgeWeight(he);
          }
        }
      }
    }

    // find edges that must be removed or inserted
    Edgelist todo_remove, todo_insert;
    for (QNodeID u = 0; u < numberOfNodes(); ++u) {
      for (const auto& delta : deltas[u]) {
        const QNodeID& v = delta.first;
        const QEdgeWeight& weight = delta.second;
        if (weight != 0) {
          if (_adjacency_matrix[u][v] == 0) {
            ASSERT(weight > 0);
            todo_insert.emplace_back(u, v);
          } else if (_adjacency_matrix[u][v] + weight == 0) {
            todo_remove.emplace_back(u, v);
          }
        }
      }
    }

    // legal movement or does it induce a cycle?
    if (!testChanges(todo_remove, todo_insert)) {
      return false; // illegal move, don't commit
    }

    // commit changes
    for (QNodeID u = 0; u < numberOfNodes(); ++u) {
      for (const auto& delta : deltas[u]) {
        const QNodeID& v = delta.first;
        const QEdgeWeight& weight = delta.second;
        _adjacency_matrix[u][v] += weight;
      }
    }
    for (const auto& edge : todo_remove) {
      const QNodeID& u = edge.first;
      const QNodeID& v = edge.second;
      _detector.disconnect(u, v);
    }
    for (const auto& edge : todo_insert) {
      const QNodeID& u = edge.first;
      const QNodeID& v = edge.second;
      _detector.connect(u, v);
    }

    return true;
  }

  /**!
   * Checks whether the quotient graph is acyclic. Note that a quotient graph that is acyclic initially should never
   * become cyclic by calling the public member functions.
   *
   * \return Whether the quotient graph is acyclic.
   */
  bool isAcyclic() const {
    std::vector<QEdgeID> rank(numberOfNodes());
    for (QNodeID u = 0; u < numberOfNodes(); ++u) {
      for (const auto& edge : _adjacency_matrix[u]) {
        const QNodeID& v = edge.first;
        const QEdgeWeight& weight = edge.second;
        if (weight != 0) {
          ASSERT(weight > 0);
          ++rank[v];
        }
      }
    }

    std::vector<QNodeID> candidates;
    for (QNodeID u = 0; u < numberOfNodes(); ++u) {
      if (rank[u] == 0) {
        candidates.push_back(u);
      }
    }

    std::size_t processed_candidates = 0;
    while (!candidates.empty()) {
      const QNodeID u = candidates.back();
      candidates.pop_back();
      ++processed_candidates;

      for (const auto& edge : _adjacency_matrix[u]) {
        const QNodeID& v = edge.first;
        const QEdgeWeight& weight = edge.second;
        if (weight != 0) {
          ASSERT(weight > 0);
          ASSERT(rank[v] > 0);
          --rank[v];
          if (rank[v] == 0) {
            candidates.push_back(v);
          }
        }
      }
    }

    return processed_candidates == numberOfNodes();
  }

  void log() const {
    LOG << "Quotient graph with" << numberOfNodes() << "nodes and" << numberOfEdges() << "directed edges";
    for (QNodeID u = 0; u < numberOfNodes(); ++u) {
      std::stringstream ss;
      ss << u << ": ";
      for (const auto& edge : _adjacency_matrix[u]) {
        const QNodeID& v = edge.first;
        const QEdgeWeight& weight = edge.second;
        ss << v << "(" << weight << "), ";
      }
      ss << "\b\b";
      LOG << ss.str();
    }
    LOG << "End of quotient graph";
  }

 private:
  using Edgelist = std::vector<std::pair<QNodeID, QNodeID>>;
  using AdjacencyMatrix = std::vector<std::unordered_map<QNodeID, QEdgeWeight>>;

  /**
   * Tests whether removing and inserting a given set of edges induces a cycle in this quotient graph.
   *
   * \param todo_remove Existing edges to be removed
   * \param todo_insert Nonexisting edges to be inserted
   * \return Whether the changes form a cycle in this quotient graph.
   */
  bool testChanges(const Edgelist& todo_remove, const Edgelist& todo_insert) {
    for (const auto& edge : todo_remove) {
      _detector.disconnect(edge.first, edge.second);
    }

    for (int i = 0; i < todo_insert.size(); ++i) {
      const auto& edge = todo_insert[i];
      if (!_detector.connect(edge.first, edge.second)) { // cycle detected, rollback
        // revert insertions
        for (int j = i - 1; j >= 0; --j) {
          const auto& rev_edge = todo_insert[j];
          _detector.disconnect(rev_edge.first, rev_edge.second);
        }

        // revert deletions
        for (const auto& rev_edge : todo_remove) {
          bool success = _detector.connect(rev_edge.first, rev_edge.second);
          ASSERT(success);
        }

        return false;
      }
    }

    return true; // no cycle
  }

  /**!
   * Construct the initial quotient graph based on a partitioned hypergraph.
   */
  void constructFromHypergraph() {
    ASSERT(_context.partition.k > 0);

    // clear previous quotient graph
    _adjacency_matrix.clear();
    _adjacency_matrix.resize(_context.partition.k);

    for (const HyperedgeID& he : _hg.edges()) {
      for (const HypernodeID& hh : _hg.heads(he)) {
        for (const HypernodeID& ht : _hg.tails(he)) {
          PartitionID from = _hg.partID(ht);
          PartitionID to = _hg.partID(hh);
          if (from != to) {
            ++_adjacency_matrix[from][to];
          }
        }
      }
    }

    // add initial edges to the cycle detector
    cycledetector::Edgelist initial_edges;
    for (QNodeID u = 0; u < numberOfNodes(); ++u) {
      for (const auto& edge : _adjacency_matrix[u]) {
        const QNodeID& v = edge.first;
        const QEdgeWeight& weight = edge.second;
        ASSERT(weight > 0);
        initial_edges.emplace_back(u, v);
      }
    }
    _detector.bulkConnect(initial_edges);
  }

  const Hypergraph& _hg;
  const Context& _context;
  AdjacencyMatrix _adjacency_matrix;
  Detector _detector;
};
} // namespace kahypar
