#pragma once

#include <vector>
#include <stack>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <set>
#include <kahypar/meta/typelist.h>

#include "kahypar/macros.h"

namespace kahypar {
namespace cycledetector {
using Edgelist = std::vector<std::pair<std::size_t, std::size_t>>;

class CycleDetector {
 public:
  /**!
   * Inserts a batch of edges without checking for cycles. Datastructure may be in an inconsistent state if
   * edges inserted using this method make the graph cyclic.
   *
   * \param edges List of edges to insert.
   */
  virtual void bulkConnect(const Edgelist& edges) = 0;

  /**!
   * Tests whether the graph remains acyclic after inserting the edge. Inserts the edge if the graph remains
   * acyclic. Does not insert the edge if it induces a cycle.
   *
   * \param s Source node, tail of the edge
   * \param t Target node, head of the edge
   * \return Whether the graph remains acyclic
   */
  virtual bool connect(std::size_t s, std::size_t t) = 0;

  /**!
   * Removes an edge from the graph.
   *
   * \param s Source node, tail of the edge
   * \param t Target node, head of the edge
   */
  virtual void disconnect(std::size_t s, std::size_t t) = 0;

  virtual void log() {}

  virtual void reset() = 0;
};

/**!
 * Base class for cycle detectors that store edges in an adjacency list.
 */
class OutCycleDetector : CycleDetector {
 public:
  void bulkConnect(const Edgelist& edges) override {
    for (const auto& pair : edges) {
      _out[pair.first].push_back(pair.second);
    }
  }

  void disconnect(std::size_t s, std::size_t t) override {
    for (std::size_t i = 0; i < _out[s].size(); ++i) {
      if (_out[s][i] == t) {
        std::swap(_out[s][i], _out[s].back());
        _out[s].pop_back();
      }
    }
  }

  void log() override {
    LOG << "State of the Cycle Detector:";
    for (std::size_t u = 0; u < _out.size(); ++u) {
      LOG << u << ":" << _out[u];
    }
    LOG << "End of the state";
  }

  void reset() override {
    const std::size_t n = _out.size();
    _out.clear();
    _out.resize(n);
  }

 protected:
  explicit OutCycleDetector(std::size_t n) :
    _out(n) {}

  /**!
   * Calculates a topological ordering of the graph.
   *
   * \return Vector vec s.t. vec[u] holds the position of node u in the topological ordering; or an empty vector
   * if the graph is no longer acyclic.
   */
  std::vector<std::size_t> calculateTopologicalOrdering() const {
    std::vector<std::size_t> rank(_out.size());
    for (const auto& neighbors : _out) {
      for (const std::size_t& v : neighbors) {
        ++rank[v];
      }
    }

    std::stack<std::size_t> candidates;
    for (std::size_t u = 0; u < _out.size(); ++u) {
      if (rank[u] == 0) {
        candidates.push(u);
      }
    }

    std::size_t pos = 0;
    std::vector<std::size_t> ordering(_out.size());
    while (!candidates.empty()) {
      std::size_t u = candidates.top();
      candidates.pop();

      for (const std::size_t& v : _out[u]) {
        ASSERT(rank[v] > 0);
        --rank[v];
        if (rank[v] == 0) {
          candidates.push(v);
        }
      }

      ordering[u] = pos++;
    }

    if (pos != _out.size()) { // not acyclic
      return {};
    }
    return ordering; // acyclic
  }

  std::vector<std::vector<std::size_t>> _out;
};

class NullCycleDetector : public CycleDetector {
 public:
  explicit NullCycleDetector(const PartitionID) {}

  void bulkConnect(const Edgelist& edges) override {
  }

  bool connect(std::size_t s, std::size_t t) override {
    return true;
  }

  void disconnect(std::size_t s, std::size_t t) override {
  }

  void reset() override {
  }
};

/**!
 * Variant that stores a topological ordering of the graph. When inserting an edge that violates the topological
 * ordering, the algorithm tests whether a new topological ordering can be obtained with the edge inserted.
 */
class KahnCycleDetector : public OutCycleDetector {
 public:
  explicit KahnCycleDetector(std::size_t n) :
    cycledetector::OutCycleDetector(n),
    _topological_ordering(n) {
    std::iota(_topological_ordering.begin(), _topological_ordering.end(), 0);
  }

  KahnCycleDetector(const KahnCycleDetector& other) = delete;
  KahnCycleDetector& operator=(const KahnCycleDetector& other) = delete;

  KahnCycleDetector(KahnCycleDetector&& other) = default;
  KahnCycleDetector& operator=(KahnCycleDetector&& other) = default;

  ~KahnCycleDetector() = default;

  void bulkConnect(const cycledetector::Edgelist& edges) override {
    OutCycleDetector::bulkConnect(edges);
    updateTopologicalOrdering();
    ASSERT(!_topological_ordering.empty());
  }

  bool connect(std::size_t s, std::size_t t) override {
    ASSERT(s < _out.size());
    ASSERT(t < _out.size());

    if (_topological_ordering[s] < _topological_ordering[t]) { // no cycle
      _out[s].push_back(t);
      return true;
    }

    _out[s].push_back(t);
    auto topological_ordering = calculateTopologicalOrdering();
    if (topological_ordering.empty()) { // cycle -> revert edge insertion
      _out[s].pop_back();
      return false;
    } else { // no cycle -> keep topological ordering and done
      _topological_ordering = topological_ordering;
      return true;
    }
  }

 private:
  std::vector<std::size_t> _topological_ordering;

  void updateTopologicalOrdering() {
    _topological_ordering = calculateTopologicalOrdering();
  }
};

/**!
 * Uses a simple DFS search to test whether new edges induce a cycle.
 */
class DFSCycleDetector : public OutCycleDetector {
 public:
  explicit DFSCycleDetector(std::size_t n) :
    OutCycleDetector(n),
    _run(n) {}

  DFSCycleDetector(const DFSCycleDetector& other) = delete;
  DFSCycleDetector& operator=(const DFSCycleDetector& other) = delete;

  DFSCycleDetector(DFSCycleDetector&& other) = default;
  DFSCycleDetector& operator=(DFSCycleDetector&& other) = default;

  ~DFSCycleDetector() = default;

  void bulkConnect(const Edgelist& edges) override {
    OutCycleDetector::bulkConnect(edges);
    //ASSERT(!calculateTopologicalOrdering().empty());
  }

  bool connect(std::size_t s, std::size_t t) override {
    ASSERT(s < _out.size());
    ASSERT(t < _out.size());

    std::stack<std::size_t> todo;
    todo.push(t);
    ++_run_id;

    while (!todo.empty()) {
      std::size_t u = todo.top();
      todo.pop();
      if (u == s) return false;
      _run[u] = _run_id;

      for (const std::size_t& v : _out[u]) {
        if (_run[v] != _run_id) {
          todo.push(v);
        }
      }
    }

    _out[s].push_back(t);
    return true;
  }

 private:
  std::vector<std::size_t> _run;
  std::size_t _run_id{0};
};

// Based on "A NEW APPROACH TO INCREMENTAL CYCLE DETECTION AND RELATED PROBLEMS" by Bender et. al.
class PseudoTopologicalOrderingCycleDetector : CycleDetector {
 public:
  explicit PseudoTopologicalOrderingCycleDetector(std::size_t n) :
    _level(n),
    _out(n),
    _in(n),
    _marked(n),
    _cbrt_sqr_order(static_cast<std::size_t>(std::cbrt(n))) {}

  PseudoTopologicalOrderingCycleDetector(const PseudoTopologicalOrderingCycleDetector& other) = delete;
  PseudoTopologicalOrderingCycleDetector& operator=(const PseudoTopologicalOrderingCycleDetector& other) = delete;

  PseudoTopologicalOrderingCycleDetector(PseudoTopologicalOrderingCycleDetector&& other) = default;
  PseudoTopologicalOrderingCycleDetector& operator=(PseudoTopologicalOrderingCycleDetector&& other) = delete;

  ~PseudoTopologicalOrderingCycleDetector() = default;

  void bulkConnect(const Edgelist& edges) override {
    _size = edges.size();
    updateDelta();
    for (const auto& edge : edges) {
      connect(edge.first, edge.second);
      --_size; // connect increments size, but we've already set size to the right value
    }
  }

  void disconnect(std::size_t s, std::size_t t) override {
    _out[s].erase(std::find(_out[s].begin(), _out[s].end(), t));
    _in[t].erase(std::find(_in[t].begin(), _in[t].end(), s));
  }

  bool connect(std::size_t u, std::size_t v) override {
    ASSERT([&](){
      for (const std::size_t& y : _out[u]) {
        if (y == v) {
          return false;
        }
      }
      return true;
    }(), "Edge already present");

    bool success = (_level[u] < _level[v])
                   ? insertEdge(u, v)
                   : performBackwardSearch(u, v);

    //LOG << "State after insertion:" << success;
    //for (std::size_t x = 0; x < _out.size(); ++x) {
    //  LOG << "Node:" << x << "Level:" << _level[x] << ":" << _out[x] << "|" << _in[x];
    //}
    ASSERT([&]() {
      for (std::size_t x = 0; x < _out.size(); ++x) {
        for (const std::size_t& y : _out[x]) {
          if (_level[x] > _level[y]) {
            return false;
          }
        }
      }
      return true;
    }(), "Pseudo topological ordering violated");
    return success;
  }

  void log() const {
    for (std::size_t x = 0; x < _out.size(); ++x) {
      LOG << "Node" << x << "| Level" << _level[x] << "| Out" << _out[x] << "| In" << _in[x];
    }
  }

  void reset() override {
    const std::size_t n = _out.size();
    _level.clear();
    _level.resize(n);
    _out.clear();
    _out.resize(n);
    _in.clear();
    _in.resize(n);
    _marked.clear();
    _marked.resize(n);
    _delta = 0;
    _size = 0;
    _mark = 0;
  }

 private:
  bool insertEdge(std::size_t u, std::size_t v) {
    _out[u].insert(v);
    if (_level[u] == _level[v]) _in[v].insert(u);
    ++_size;
    updateDelta();
    return true;
  }

  void updateDelta() {
    _delta = std::min(static_cast<std::size_t>(std::sqrt(_size)), _cbrt_sqr_order);
  }

  bool performBackwardSearch(std::size_t u, std::size_t v) {
    ASSERT(_level[u] >= _level[v]);

    std::size_t visited = 0;
    ++_mark;
    std::stack<std::size_t> todo({u}); // LIFO for depth first search

    while (!todo.empty()) {
      const std::size_t y = todo.top();
      todo.pop();

      if (_marked[y] != _mark) {
        ++visited;
        _marked[y] = _mark;

        // cycle detected
        if (y == v) {
          return false;
        }

        // add unmarked predecessors to the to-do LIFO
        for (const std::size_t& x : _in[y]) {
          if (_marked[x] != _mark && _level[x] == _level[y]) {
            todo.push(x);
          }
        }

        // cancel backward search after at least _delta nodes were visited
        if (visited >= _delta) {
          break;
        }
      }
    }

    // easy case -- structure remains in good state
    if (visited < _delta && _level[u] == _level[v]) {
      return insertEdge(u, v);
    }

    // remaining two cases -- there might be a cycle, ordering must be updated
    ASSERT(visited >= _delta || _level[u] > _level[v]);
    std::size_t old_level_v = _level[v];
    _level[v] = (visited < _delta) ? _level[u] : _level[u] + 1;

    if (visited >= _delta) { // backward search was inconclusive, reset B (_marked in code)
      for (std::size_t& mark : _marked) {
        mark = 0;
      }
      _marked[u] = _mark;
    }

    // we split Step 3 into two steps: the first one does not modify the _in sets but only searches for a cycle, the
    // second part updates the _in sets -- this way, we don't have to rollback changes made to the _in sets
    if (!performForwardCycleSearch(u, v)) {
      _level[v] = old_level_v;
      return false;
    }

    // we're in the clear now, edge insertion is ok, no more rollbacks
    _in[v].clear();
    performForwardSearch(u, v);
    return true;
  }

  bool performForwardCycleSearch(std::size_t u, std::size_t v) {
    std::stack<std::pair<std::size_t, std::size_t>> todo;
    for (const std::size_t& y : _out[v]) {
      todo.push({v, y});
    }

    std::vector<std::pair<std::size_t, std::size_t>> changed_levels;
    bool success = true;

    while (!todo.empty()) {
      const auto edge = todo.top();
      todo.pop();
      const std::size_t x = edge.first;
      const std::size_t y = edge.second;

      if (_marked[y] == _mark) { // found cycle
        success = false;
        break;
      }

      if (_level[x] > _level[y]) {
        changed_levels.emplace_back(y, _level[y]);
        _level[y] = _level[x];

        for (const std::size_t& y_prime : _out[y]) {
          todo.push({y, y_prime});
        }
      }
    }

    // rollback
    for (const auto& changed_level : changed_levels) {
      const std::size_t& x = changed_level.first;
      const std::size_t& level = changed_level.second;
      _level[x] = level;
    }

    return success;
  }

  void performForwardSearch(std::size_t u, std::size_t v) {
    std::stack<std::pair<std::size_t, std::size_t>> todo;
    for (const std::size_t& y : _out[v]) {
      todo.push({v, y});
    }

    while (!todo.empty()) {
      const auto edge = todo.top();
      todo.pop();
      const std::size_t x = edge.first;
      const std::size_t y = edge.second;

      if (_level[x] == _level[y]) {
        _in[y].insert(x);
      } else if (_level[x] > _level[y]) {
        _level[y] = _level[x];
        _in[y].clear();
        _in[y].insert(x);
        for (const std::size_t& y_prime : _out[y]) {
          todo.push({y, y_prime});
        }
      }
    }

    insertEdge(u, v);
  }

  std::vector<std::size_t> _level;
  std::vector<std::set<std::size_t>> _out;
  std::vector<std::set<std::size_t>> _in;
  std::vector<std::size_t> _marked;
  const std::size_t _cbrt_sqr_order;

  std::size_t _delta{0};
  std::size_t _size{0};
  std::size_t _mark{0};
};
} // namespace cycledetector

// lift public interface into kahypar namespace
using cycledetector::CycleDetector;
using cycledetector::KahnCycleDetector;
using cycledetector::DFSCycleDetector;
using cycledetector::PseudoTopologicalOrderingCycleDetector;
using cycledetector::NullCycleDetector;

using CycleDetectors = meta::Typelist<KahnCycleDetector, DFSCycleDetector, PseudoTopologicalOrderingCycleDetector>;
} // namespace kahypar