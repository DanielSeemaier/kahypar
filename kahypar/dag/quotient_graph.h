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
namespace ds {
class DynamicUnweightedGraph {
 public:
  using AdjacencyList = std::vector<std::vector<QNodeID>>;

  explicit DynamicUnweightedGraph(QNodeID numberOfNodes) :
    _outs(numberOfNodes),
    _ins(numberOfNodes) {}

  DynamicUnweightedGraph(AdjacencyList outs, AdjacencyList ins) :
    _outs(std::move(outs)),
    _ins(std::move(ins)) {
    ASSERT(_outs.size() == _ins.size());
  }

  QNodeID numberOfNodes() const {
    return _outs.size();
  }

  QNodeID newNode() {
    _ins.emplace_back();
    _outs.emplace_back();
    return numberOfNodes() - 1;
  }

  void newDirectedEdge(QNodeID from, QNodeID to) {
    ASSERT(from < numberOfNodes());
    ASSERT(to < numberOfNodes());
    ASSERT(std::find(_outs[from].begin(), _outs[from].end(), to) == _outs[from].end());
    ASSERT(std::find(_ins[to].begin(), _ins[to].end(), from) == _ins[to].end());
    _outs[from].push_back(to);
    _ins[to].push_back(from);
  }

  void newUndirectedEdge(QNodeID u, QNodeID v) {
    newDirectedEdge(u, v);
    newDirectedEdge(v, u);
  }

  QNodeID outDegree(QNodeID node) const {
    return outs(node).size();
  }

  const std::vector<QNodeID>& outs(QNodeID node) const {
    ASSERT(node < _outs.size());
    return _outs[node];
  }

  QNodeID inDegree(QNodeID node) const {
    return ins(node).size();
  }

  const std::vector<QNodeID>& ins(QNodeID node) const {
    ASSERT(node < _ins.size());
    return _ins[node];
  }

  void log() const {
    LOG << "Dynamic unweighted graph with" << numberOfNodes() << "nodes";
    for (QNodeID u = 0; u < numberOfNodes(); ++u) {
      LOG << "[" << u << "] outs =" << _outs[u] << ", ins =" << _ins[u];
    }
    LOG << "End of dyanmic unweighted graph";
  }

 private:
  AdjacencyList _outs;
  AdjacencyList _ins;
};

using QuotientMoveGraph = DynamicUnweightedGraph;

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
    bool changed;
    return update(hn, from, to, changed);
  }

  bool update(const HypernodeID hn, const PartitionID from, const PartitionID to, bool& changed) {
    if (from == to) {
      changed = false;
      return true;
    }

    AdjacencyMatrix deltas(numberOfNodes());

    // find changes to edge weights of the quotient graph caused by this movement
    for (const HyperedgeID& he : _hg.incidentEdges(hn)) {
      if (_hg.isHead(hn, he)) { // hn is head
        for (const HypernodeID& ht : _hg.tails(he)) {
          if (_hg.partID(ht) != from) {
            deltas[_hg.partID(ht)][from] -= 1;
          }
          if (_hg.partID(ht) != to) {
            deltas[_hg.partID(ht)][to] += 1;
          }
        }
      } else { // hn is tail
        ASSERT(_hg.isTail(hn, he));
        for (const HypernodeID& hh : _hg.heads(he)) {
          if (_hg.partID(hh) != from) {
            deltas[from][_hg.partID(hh)] -= 1;
          }
          if (_hg.partID(hh) != to) {
            deltas[to][_hg.partID(hh)] += 1;
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
    // note: if successful, testChanges does not rollback edge insertions / deletions, i.e. the detector is up to date
    // with the changes after a successful call
    if (!testChanges(todo_remove, todo_insert)) {
      changed = false;
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

    changed = !todo_insert.empty() || !todo_remove.empty();
    return true;
  }

  std::vector<QNodeID> computeWeakTopologicalOrdering() const {
    return computeTopologicalOrdering(false);
  }

  std::vector<QNodeID> computeStrictTopologicalOrdering() const {
    return computeTopologicalOrdering(true);
  }

  std::vector<QNodeID> computeTopologicalOrdering(bool strict = true) const {
    std::vector<QEdgeID> in_degree(numberOfNodes());
    for (QNodeID u = 0; u < numberOfNodes(); ++u) {
      for (const auto& edge : _adjacency_matrix[u]) {
        const QNodeID& v = edge.first;
        const QEdgeWeight& weight = edge.second;
        if (weight != 0) {
          ASSERT(weight > 0, V(u) << V(v) << V(weight));
          ++in_degree[v];
        }
      }
    }

    std::vector<QNodeID> todo;
    for (QNodeID u = 0; u < numberOfNodes(); ++u) {
      if (in_degree[u] == 0) {
        todo.push_back(u);
      }
    }

    std::vector<QNodeID> weak_topological_ordering(numberOfNodes());
    QNodeID current_level = 0;
    std::size_t processed_nodes = 0;
    while (!todo.empty()) {
      std::vector<QNodeID> next_todo;
      for (const QNodeID& u : todo) {
        weak_topological_ordering[u] = current_level;
        if (strict) {
          ++current_level;
        }
        ++processed_nodes;

        for (const auto& edge : _adjacency_matrix[u]) {
          const QNodeID& v = edge.first;
          const QEdgeWeight& weight = edge.second;
          if (weight != 0) {
            ASSERT(weight > 0);
            ASSERT(in_degree[v] > 0);
            --in_degree[v];
            if (in_degree[v] == 0) {
              next_todo.push_back(v);
            }
          }
        }
      }

      if (!strict) {
        ++current_level;
      }
      todo = next_todo;
    }

    if (processed_nodes != numberOfNodes()) {
      return {};
    } else {
      return weak_topological_ordering;
    }
  }

  QuotientMoveGraph computeMoveGraph() const {
    QuotientMoveGraph qmg(numberOfNodes());
    const auto ordering = computeWeakTopologicalOrdering();

    for (QNodeID u = 0; u < numberOfNodes(); ++u) {
      std::vector<QNodeID> min_neighbors;
      for (const auto& edge : _adjacency_matrix[u]) {
        const QNodeID v = edge.first;
        const QEdgeWeight& weight = edge.second;
        if (weight == 0) {
          continue;
        }

        if (min_neighbors.empty()) {
          min_neighbors.push_back(v);
        } else {
          const QNodeID min_neighbor = min_neighbors[0];
          if (ordering[v] == ordering[min_neighbor]) {
            min_neighbors.push_back(v);
          } else if (ordering[v] < ordering[min_neighbor]) {
            min_neighbors.clear();
            min_neighbors.push_back(v);
          }
        }
      }

      for (const QNodeID& min_neighbor : min_neighbors) {
        if (ordering[min_neighbor] == ordering[u] + 1) {
          qmg.newUndirectedEdge(u, min_neighbor);
        } else {
          qmg.newDirectedEdge(u, min_neighbor);
        }
      }
    }

    return qmg;
  }

  /**!
   * Checks whether the quotient graph is acyclic. Note that a quotient graph that is acyclic initially should never
   * become cyclic by calling the public member functions, i.e. this is just a method for ASSERTs.
   *
   * \return Whether the quotient graph is acyclic.
   */
  bool isAcyclic() const {
    return computeTopologicalOrdering().size() == numberOfNodes();
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
} // namespace ds

using ds::QuotientGraph;
using ds::QuotientMoveGraph;
} // namespace kahypar
