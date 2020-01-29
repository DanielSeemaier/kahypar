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
#include <queue>

#include <boost/iterator/filter_iterator.hpp>

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

  bool findPath(QNodeID from, QNodeID to, std::vector<QNodeID>& path) const {
    ASSERT(from < numberOfNodes());
    ASSERT(to < numberOfNodes());
    ASSERT(path.empty());

    std::vector<QNodeID> parent(numberOfNodes());
    std::iota(parent.begin(), parent.end(), 0);

    std::queue<QNodeID> todo;
    todo.push(from);
    while (!todo.empty()) {
      QNodeID u = todo.front();
      todo.pop();

      for (const QNodeID& v  : outs(u)) {
        if (v != from && parent[v] == v) {
          parent[v] = u;
          if (v == to) {
            goto done;
          } else {
            todo.push(v);
          }
        }
      }
    }
    done:

    if (parent[to] != to) {
      for (QNodeID v = to; v != from; v = parent[v]) {
        path.push_back(v);
      }
      path.push_back(from);
      std::reverse(path.begin(), path.end());
      return true;
    }
    return false;
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

class QuotientGraph;

struct IsNonzeroEdge {
  explicit IsNonzeroEdge(const QuotientGraph& qg, const QNodeID u) :
    _qg(qg),
    _u(u) {}

  bool operator()(QNodeID v);

 private:
  const QuotientGraph& _qg;
  const QNodeID _u;
};

using QuotientGraphOutsIterator = boost::filter_iterator<IsNonzeroEdge, const QNodeID*>;

struct OutsIteratorPair {
  QuotientGraphOutsIterator begin() {
    return _begin;
  }

  QuotientGraphOutsIterator end() {
    return _end;
  }

  QuotientGraphOutsIterator _begin;
  QuotientGraphOutsIterator _end;
};

class QuotientGraph {
 public:
  explicit QuotientGraph(const PartitionID k) :
    _nodes(k) {
    std::iota(_nodes.begin(), _nodes.end(), 0);
  }

  virtual ~QuotientGraph() = default;

  virtual QNodeID numberOfNodes() const = 0;
  virtual QEdgeID numberOfEdges() const = 0;
  virtual QEdgeWeight edgeWeight(QNodeID u, QNodeID v) const = 0;
  virtual bool testAndUpdateBeforeMovement(HypernodeID hn, PartitionID from_part, PartitionID to_part) = 0;
  virtual void updateAfterMovement(HypernodeID hn, PartitionID from_part, PartitionID to_part) = 0;

  virtual void updateBeforeMovement(HypernodeID hn, PartitionID from_part, PartitionID to_part) {
    updateAfterMovement(hn, from_part, to_part);
  }

  virtual bool connected(QNodeID u, QNodeID v) const = 0;
  virtual OutsIteratorPair outs(QNodeID u) const = 0;
  virtual void preUncontraction(HypernodeID representant) = 0;
  virtual void postUncontraction(HypernodeID representant, const std::vector<HypernodeID>&& partners) = 0;

  virtual const std::vector<QNodeID>& topologicalOrdering() const {
    if (!cachedTopologicalOrderingStillTopological()) {
      auto* nc_this = const_cast<QuotientGraph*>(this);
      nc_this->_cached_topological_ordering = computeTopologicalOrdering();
      if (!_cached_topological_ordering.empty()) {
        nc_this->_cached_inverse_topological_ordering.resize(numberOfNodes());
        for (QNodeID u = 0; u < numberOfNodes(); ++u) {
          nc_this->_cached_inverse_topological_ordering[_cached_topological_ordering[u]] = u;
        }
      }
      nc_this->_outdated_topological_ordering = false;
      nc_this->_changed = true;
    }
    return _cached_topological_ordering;
  }

  virtual const std::vector<QNodeID>& inverseTopologicalOrdering() const {
    if (!cachedTopologicalOrderingStillTopological()) {
      topologicalOrdering();
    }
    return _cached_inverse_topological_ordering;
  }

  virtual std::vector<QNodeID> randomizedTopologicalOrdering() const {
    return computeTopologicalOrdering(true, true);
  }

  virtual const std::vector<QNodeID>& weakTopologicalOrdering() const {
    if (!cachedWeakTopologicalOrderingStillWeakTopological()) {
      const_cast<QuotientGraph*>(this)->_cached_weak_topological_ordering = computeTopologicalOrdering(false);
      const_cast<QuotientGraph*>(this)->_outdated_weak_topological_ordering = false;
      const_cast<QuotientGraph*>(this)->_changed = true;
    }
    return _cached_weak_topological_ordering;
  }

  virtual void log() const {
    LOG << "Quotient graph with" << numberOfNodes() << "nodes and" << numberOfEdges() << "directed edges";
    for (QNodeID u = 0; u < numberOfNodes(); ++u) {
      std::stringstream ss;
      ss << u << ": ";
      for (const auto& v : outs(u)) {
        ss << v << ", ";
      }
      ss << "\b\b";
      LOG << ss.str();
    }
    LOG << "End of quotient graph";
  }

  virtual bool isAcyclic() const {
    return !topologicalOrdering().empty();
  }

  virtual void rebuild() = 0;

  virtual bool changed() {
    return _changed;
  }

  virtual void resetChanged() {
    _changed = false;
  }

 protected:
  void notifyGraphChanged() {
    _outdated_topological_ordering = true;
    _outdated_weak_topological_ordering = true;
  }

  std::vector<QNodeID> computeTopologicalOrdering(bool strict = true, bool randomized = false) const {
    ASSERT(!randomized, "randomized topological ordering currently not implemented");
    std::vector<QEdgeID> in_degree(numberOfNodes());
    for (QNodeID u = 0; u < numberOfNodes(); ++u) {
      for (const QNodeID& v : outs(u)) {
        ++in_degree[v];
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

        for (const QNodeID& v : outs(u)) {
          ASSERT(in_degree[v] > 0);
          --in_degree[v];
          if (in_degree[v] == 0) {
            next_todo.push_back(v);
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

  bool cachedTopologicalOrderingStillTopological() const {
    if (_cached_topological_ordering.empty()) {
      return false;
    }
    if (!_outdated_topological_ordering) {
      return true;
    }

    for (QNodeID u = 0; u < numberOfNodes(); ++u) {
      for (const QNodeID v : outs(u)) {
        if (_cached_topological_ordering[u] >= _cached_topological_ordering[v]) {
          return false;
        }
      }
    }

    auto* nc_this = const_cast<QuotientGraph*>(this);
    nc_this->_outdated_topological_ordering = false;
    return true;
  }

  bool cachedWeakTopologicalOrderingStillWeakTopological() const {
    if (_cached_weak_topological_ordering.empty()) {
      return false;
    }
    if (!_outdated_weak_topological_ordering) {
      return true;
    }

    for (QNodeID u = 0; u < numberOfNodes(); ++u) {
      for (const QNodeID v : outs(u)) {
        if (_cached_weak_topological_ordering[u] > _cached_weak_topological_ordering[v]) {
          return false;
        }
      }
    }

    auto* nc_this = const_cast<QuotientGraph*>(this);
    nc_this->_outdated_weak_topological_ordering = false;
    return true;
  }

 protected:
  std::vector<QNodeID> _nodes;

 private:
  bool _outdated_topological_ordering{false};
  bool _outdated_weak_topological_ordering{false};
  bool _changed{false};
  std::vector<QNodeID> _cached_topological_ordering{};
  std::vector<QNodeID> _cached_weak_topological_ordering{};
  std::vector<QNodeID> _cached_inverse_topological_ordering{};
};

bool IsNonzeroEdge::operator()(QNodeID v) {
  return _qg.connected(_u, v);
}

template<typename CycleDetector>
class AdjacencyMatrixQuotientGraph : public QuotientGraph {
  using AdjacencyMatrix = std::vector<std::vector<QEdgeWeight>>;
  using Edgelist = std::vector<std::pair<QNodeID, QNodeID>>;

 public:
  AdjacencyMatrixQuotientGraph(const Hypergraph& hg, const Context& context) :
    QuotientGraph(context.partition.k),
    _hg(hg),
    _context(context),
    _detector(context.partition.k),
    _adj_matrix(context.partition.k, std::vector<QEdgeWeight>(context.partition.k)),
    _delta_matrix(context.partition.k, std::vector<QEdgeWeight>(context.partition.k)),
    _forbidden_edges_cache(context.partition.k, std::vector<bool>(context.partition.k)) {
    buildFromHypergraph();
  }

  AdjacencyMatrixQuotientGraph(const AdjacencyMatrixQuotientGraph&) = delete;
  AdjacencyMatrixQuotientGraph& operator= (const AdjacencyMatrixQuotientGraph&) = delete;

  AdjacencyMatrixQuotientGraph(AdjacencyMatrixQuotientGraph&&) = delete;
  AdjacencyMatrixQuotientGraph& operator= (AdjacencyMatrixQuotientGraph&&) = delete;

  ~AdjacencyMatrixQuotientGraph() override = default;

  QNodeID numberOfNodes() const override {
    return _context.partition.k;
  }

  QEdgeID numberOfEdges() const override {
    return _num_edges;
  }

  QEdgeWeight edgeWeight(const QNodeID u, const QNodeID v) const override {
    return _adj_matrix[u][v];
  }

  bool connected(const QNodeID u, const QNodeID v) const override {
    return edgeWeight(u, v) > 0;
  }

  bool testAndUpdateBeforeMovement(const HypernodeID hn, const PartitionID from_part, const PartitionID to_part) override {
    ASSERT(_hg.partID(hn) == from_part, "You must call this method before moving the hypernode!");
    ASSERT(from_part != to_part);

//    if (_forbidden_edges_cache[from_part][to_part]) {
//      return false;
//    }

    resetDeltaMatrix();

    // note: during graph construction, edges are counted twice, thus we increment / decrement the delta matrix by 2

    // changes due to tail --> hn edges
    for (const HyperedgeID& he : _hg.incidentHeadEdges(hn)) {
      for (const HypernodeID& tail : _hg.tails(he)) {
        if (_hg.partID(tail) != from_part) {
          changeDeltaMatrixBy(_hg.partID(tail), from_part, -2);
        }
        if (_hg.partID(tail) != to_part) {
          changeDeltaMatrixBy(_hg.partID(tail), to_part, 2);
        }
      }
    }

    // changes to to hn --> head edges
    for (const HyperedgeID& he : _hg.incidentTailEdges(hn)) {
      for (const HypernodeID& head : _hg.heads(he)) {
        if (_hg.partID(head) != from_part) {
          changeDeltaMatrixBy(from_part, _hg.partID(head), -2);
        }
        if (_hg.partID(head) != to_part) {
          changeDeltaMatrixBy(to_part, _hg.partID(head), 2);
        }
      }
    }

    // find edges that must be removed or inserted
    Edgelist todo_remove, todo_insert;
    for (const auto& entry : _used_delta_entries) {
      const QNodeID u = entry.first;
      const QNodeID v = entry.second;
      ASSERT(_adj_matrix[u][v] >= 0);
      ASSERT(_adj_matrix[u][v] + _delta_matrix[u][v] >= 0);
      if (_adj_matrix[u][v] > 0 && _adj_matrix[u][v] + _delta_matrix[u][v] == 0) {
        todo_remove.emplace_back(u, v);
      } else if (_adj_matrix[u][v] == 0 && _delta_matrix[u][v] > 0) {
        todo_insert.emplace_back(u, v);
      }
    }

    // legal movement or does it induce a cycle?
    // note: if successful, testChanges does not rollback edge insertions / deletions, i.e. the detector is up to date
    // with the changes after a successful call
    if (!testChanges(todo_remove, todo_insert)) {
      cacheForbiddenEdge(from_part, to_part);
      return false; // illegal move, don't commit
    }

    clearForbiddenEdgesCache();
    commitDeltaMatrix();

     if (!todo_insert.empty()) {
       notifyGraphChanged();
     }
    return true;
  }

  void updateAfterMovement(const HypernodeID hn, const PartitionID from_part, const PartitionID to_part) override {
    ASSERT(false, "Not implemented");
  }

  OutsIteratorPair outs(const QNodeID u) const override {
    IsNonzeroEdge predicate(*this, u);
    QuotientGraphOutsIterator first(predicate, _nodes.data(), _nodes.data() + numberOfNodes());
    QuotientGraphOutsIterator last(predicate, _nodes.data() + numberOfNodes(), _nodes.data() + numberOfNodes());
    return {first, last};
  }

  void preUncontraction(const HypernodeID representant) override {
    removeHypernodePreUncontraction(representant);
  }

  void postUncontraction(const HypernodeID representant, const std::vector<HypernodeID>&& partners) override {
    ASSERT(partners.size() == 1);
    const HypernodeID partner = partners.front();

    addHypernodePostUncontraction(representant, partner);
    addHypernodePostUncontraction(partner, representant);

    ASSERT([&]() {
      ASSERT_THAT_ADJACENCY_MATRIX_IS_CORRECT();
      return true;
    }());

    clearForbiddenEdgesCache();
  }

  void rebuild() override {
    _adj_matrix.clear();
    _adj_matrix.resize(_context.partition.k, std::vector<QEdgeWeight>(_context.partition.k));
    _num_edges = 0;
    _detector.reset();
    buildFromHypergraph();
  }

 private:
  void buildFromHypergraph() {
    for (const HypernodeID& hn : _hg.nodes()) {
      addHypernode(hn);
    }

    ASSERT([&]() {
      ASSERT_THAT_ADJACENCY_MATRIX_IS_CORRECT();
      return true;
    }());

    clearForbiddenEdgesCache();
  }

  void removeHypernodePreUncontraction(const HypernodeID hn) {
    std::vector<std::pair<QNodeID, QNodeID>> edges;

    // hn --> head
    for (const HyperedgeID& he : _hg.incidentTailEdges(hn)) {
      for (const HypernodeID& head : _hg.heads(he)) {
        const PartitionID tail_part = _hg.partID(hn);
        const PartitionID head_part = _hg.partID(head);
        if (tail_part != head_part) {
          ASSERT(_adj_matrix[tail_part][head_part] >= 2);
          if (_adj_matrix[tail_part][head_part] == 2) {
            edges.emplace_back(tail_part, head_part);
          }
          _adj_matrix[tail_part][head_part] -= 2;
        }
      }
    }

    // tail --> hn
    for (const HyperedgeID& he : _hg.incidentHeadEdges(hn)) {
      for (const HypernodeID& tail : _hg.tails(he)) {
        const PartitionID tail_part = _hg.partID(tail);
        const PartitionID head_part = _hg.partID(hn);
        if (tail_part != head_part) {
          ASSERT(_adj_matrix[tail_part][head_part] >= 2);
          if (_adj_matrix[tail_part][head_part] == 2) {
            edges.emplace_back(tail_part, head_part);
          }
          _adj_matrix[tail_part][head_part] -= 2;
        }
      }
    }

    for (const auto& edge : edges) {
      _detector.disconnect(edge.first, edge.second);
    }
  }

  void addHypernode(const HypernodeID hn) {
    std::vector<std::pair<std::size_t, std::size_t>> edges;

    // hn --> head
    for (const HyperedgeID& he : _hg.incidentTailEdges(hn)) {
      for (const HypernodeID& head : _hg.heads(he)) {
        const PartitionID tail_part = _hg.partID(hn);
        const PartitionID head_part = _hg.partID(head);
        if (tail_part != head_part) {
          if (_adj_matrix[tail_part][head_part] == 0) {
            edges.emplace_back(tail_part, head_part);
          }
          ++_adj_matrix[tail_part][head_part];
        }
      }
    }

    // tail --> hn
    for (const HyperedgeID& he : _hg.incidentHeadEdges(hn)) {
      for (const HypernodeID& tail : _hg.tails(he)) {
        const PartitionID tail_part = _hg.partID(tail);
        const PartitionID head_part = _hg.partID(hn);
        if (tail_part != head_part) {
          if (_adj_matrix[tail_part][head_part] == 0) {
            edges.emplace_back(tail_part, head_part);
          }
          ++_adj_matrix[tail_part][head_part];
        }
      }
    }

    if (!edges.empty()) {
      notifyGraphChanged();
      _detector.bulkConnect(edges);
      _num_edges += edges.size();
    }
  }

  void addHypernodePostUncontraction(const HypernodeID hn, const HypernodeID partner) {
    std::vector<std::pair<std::size_t, std::size_t>> edges;

    // hn --> head
    for (const HyperedgeID& he : _hg.incidentTailEdges(hn)) {
      for (const HypernodeID& head : _hg.heads(he)) {
        const PartitionID tail_part = _hg.partID(hn);
        const PartitionID head_part = _hg.partID(head);
        if (tail_part != head_part) {
          if (_adj_matrix[tail_part][head_part] == 0) {
            edges.emplace_back(tail_part, head_part);
          }
          _adj_matrix[tail_part][head_part] += (head == partner) ? 1 : 2;
        }
      }
    }

    // tail --> hn
    for (const HyperedgeID& he : _hg.incidentHeadEdges(hn)) {
      for (const HypernodeID& tail : _hg.tails(he)) {
        const PartitionID tail_part = _hg.partID(tail);
        const PartitionID head_part = _hg.partID(hn);
        if (tail_part != head_part) {
          if (_adj_matrix[tail_part][head_part] == 0) {
            edges.emplace_back(tail_part, head_part);
          }
          _adj_matrix[tail_part][head_part] += (tail == partner) ? 1 : 2;
        }
      }
    }

    if (!edges.empty()) {
      notifyGraphChanged();
      _detector.bulkConnect(edges);
      _num_edges += edges.size();
    }
  }

  void changeDeltaMatrixBy(const std::size_t u, const std::size_t v, const QEdgeWeight delta) {
    if (_delta_matrix[u][v] == 0) {
      _used_delta_entries.emplace_back(u, v);
    }
    _delta_matrix[u][v] += delta;
  }

  void resetDeltaMatrix() {
    for (const auto& entry : _used_delta_entries) {
      _delta_matrix[entry.first][entry.second] = 0;
    }
    _used_delta_entries.clear();
  }

  void commitDeltaMatrix() {
    for (const auto& entry : _used_delta_entries) {
      const QNodeID u = entry.first;
      const QNodeID v = entry.second;
      _adj_matrix[u][v] += _delta_matrix[u][v];
      _delta_matrix[u][v] = 0;
    }
  }

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

  void cacheForbiddenEdge(std::size_t u, std::size_t v) {
    if (!_forbidden_edges_cache[u][v]) {
      _forbidden_edges_cache[u][v] = true;
      _forbidden_edges_cache_entries.emplace_back(u, v);
    }
  }

  void clearForbiddenEdgesCache() {
    for (const auto& entry : _forbidden_edges_cache_entries) {
      _forbidden_edges_cache[entry.first][entry.second] = false;
    }
    _forbidden_edges_cache_entries.clear();
  }

  void ASSERT_THAT_ADJACENCY_MATRIX_IS_CORRECT() {
    AdjacencyMatrix tmp_matrix(numberOfNodes(), std::vector<QEdgeWeight>(numberOfNodes()));

    for (const HypernodeID& hn : _hg.nodes()) {
      for (const HyperedgeID& he : _hg.incidentTailEdges(hn)) {
        for (const HypernodeID& head : _hg.heads(he)) {
          const PartitionID tail_part = _hg.partID(hn);
          const PartitionID head_part = _hg.partID(head);
          if (tail_part != head_part) {
            ++tmp_matrix[tail_part][head_part];
          }
        }
      }

      // tail --> hn
      for (const HyperedgeID& he : _hg.incidentHeadEdges(hn)) {
        for (const HypernodeID& tail : _hg.tails(he)) {
          const PartitionID tail_part = _hg.partID(tail);
          const PartitionID head_part = _hg.partID(hn);
          if (tail_part != head_part) {
            ++tmp_matrix[tail_part][head_part];
          }
        }
      }
    }

    for (QNodeID u = 0; u < numberOfNodes(); ++u) {
      for (QNodeID v = 0; v < numberOfNodes(); ++v) {
        ASSERT(tmp_matrix[u][v] == _adj_matrix[u][v], V(tmp_matrix[u][v]) << V(_adj_matrix[u][v]) << V(u) << V(v));
      }
    }
  }

  const Hypergraph& _hg;
  const Context& _context;
  CycleDetector _detector;
  std::vector<std::vector<bool>> _forbidden_edges_cache;
  std::vector<std::pair<std::size_t, std::size_t>> _forbidden_edges_cache_entries;

  AdjacencyMatrix _adj_matrix;
  AdjacencyMatrix _delta_matrix;
  std::vector<std::pair<std::size_t, std::size_t>> _used_delta_entries{};
  QEdgeID _num_edges{0};
};
} // namespace ds

using ds::QuotientGraph;
using ds::AdjacencyMatrixQuotientGraph;
using CyclicQuotientGraph = AdjacencyMatrixQuotientGraph<NullCycleDetector>;
} // namespace kahypar
