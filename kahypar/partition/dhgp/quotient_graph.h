#pragma once

#ifndef KAHYPAR_ENABLE_DHGP
#error This file should not be included when building without KAHYPAR_ENABLE_DHGP
#endif// KAHYPAR_ENABLE_DHGP

namespace kahypar::dhgp {
template<typename CycleDetector>
class QuotientGraph {
 public:
  using QNodeID = PartitionID;

  struct QEdge {
    QNodeID u, v;
  };

  QuotientGraph(const Hypergraph& hg, const Context& context)
      : _hg(hg),
        _context(context),
        _cycle_detector(_hg, _context) {}

  void initialize() {
    rebuildAdjacencyMatrix();
    rebuildHNConnectivity();
    flushBlockedMoves();
    recomputeTopologicalOrdering();
  }

  bool testAndUpdateBeforeMovement(const HypernodeID hn, const PartitionID to) {
    if (blockedMove(hn, to)) return false;
    return false;
  }

  const std::vector<QNodeID>& topologicalOrdering() const {
    if (_topord_outdated) {
      const_cast<QuotientGraph*>(this)->recomputeTopologicalOrdering();
    }
    ASSERT(!_topord_outdated);
    return _cached_topord;
  }

  const std::vector<QNodeID>& invertedTopologicalOrdering() const {
    if (_topord_outdated) {
      const_cast<QuotientGraph*>(this)->recomputeTopologicalOrdering();
    }
    ASSERT(!_topord_outdated);
    return _cached_inverted_topord;
  }

 private:
  void rebuildAdjacencyMatrix() {
    const PartitionID k = _context.partition.k;
    _adj_matrix.clear();
    _adj_matrix.resize(k * k);

    for (const HypernodeID& hn : _hg.nodes()) {
      const PartitionID part_hn = _hg.partID(hn);

      for (const HyperedgeID& he : _hg.incidentHeadEdges(hn)) {
        for (PartitionID part = 0; part < k; ++part) {
          const HypernodeID min_pins_in_part = (part == part_hn) ? 2 : 1;
          if (_hg.pinCountInPart(he, part) >= min_pins_in_part) {
            ++adjMatrix(part, part_hn);
          }
        }
      }
    }
  }

  void rebuildHNConnectivity() {
    const PartitionID k = _context.partition.k;
    const HypernodeID n = _hg.initialNumNodes();
    _hn_connectivity.clear();
    _hn_connectivity.resize(n * k);

    for (const HyperedgeID& he : _hg.edges()) {
      for (const HypernodeID& u : _hg.tailPins(he)) {
        for (const HypernodeID& v : _hg.headPins(he)) {
          const PartitionID part_u = _hg.partID(u);
          const PartitionID part_v = _hg.partID(v);

          // u --> v is cut edge
          ++hnConnectivityTo(u, part_v);
          ++hnConnectivityTo(v, part_u);
        }
      }
    }
  }

  std::vector<QEdge> determineEdgeRemovalsDueToMovement(const HypernodeID hn, const PartitionID to) const {
    const PartitionID from = _hg.partID(hn);
    const auto& topord = topologicalOrdering();
    const auto& inverted_topord = invertedTopologicalOrdering();

    std::vector<QEdge> edges;

    for (QNodeID i = 0; i < inverted_topord[from]; ++i) {
      const QNodeID u = topord[i]; // u --> from is potentially a QG edge

      ASSERT(from != u);
      ASSERT(inverted_topord[u] < inverted_topord[from]);

      if (hnConnectivityTo(from, u) > 0 // there are tail --> hn edges in part u
          && adjMatrix(u, from) == 1) { // and exactly on head in from induces the QG edge
        edges.emplace_back(u, from);
      }
    }

    return edges;
  }

  void determineEdgeInsertionsDueToMovement(const HypernodeID hn, const PartitionID to) const {

  }

  HypernodeID& hnConnectivityTo(const HypernodeID hn, const PartitionID part) {
    return _hn_connectivity[hn * _context.partition.k + part];
  }

  HypernodeID& adjMatrix(const HypernodeID u, const HypernodeID v) {
    return _adj_matrix[u * _context.partition.k + v];
  }

  auto& blockedMove(const HypernodeID hn, const PartitionID part) {
    return _blocked_moves[hn * _context.partition.k + part];
  }

  void flushBlockedMoves() {
    _blocked_moves.clear();
    _blocked_moves.resize(_hg.initialNumNodes() * _context.partition.k);
  }

  void invalidateCachedTopologicalOrdering() {
    _topord_outdated = true;
  }

  void recomputeTopologicalOrdering() {
    const PartitionID k = _context.partition.k;

    std::vector<QNodeID> indegree(k);
    for (QNodeID u = 0; u < k; ++u) {
      for (QNodeID v = 0; v < k; ++v) {
        if (u == v) continue; // self loop
        if (adjMatrix(u, v) == 0) continue; // zero weight edge
        ++indegree[v];
      }
    }

    std::vector<QNodeID> candidates;
    for (QNodeID u = 0; u < k; ++u) {
      if (indegree[u] > 0) continue; // not a root
      candidates.push_back(u);
    }

    _cached_topord.clear();
    _cached_topord.reserve(k);
    while (!candidates.empty()) {
      const QNodeID u = candidates.front();
      std::swap(candidates.front(), candidates.back());
      candidates.pop_back();

      _cached_topord.push_back(u);

      for (QNodeID v = 0; v < k; ++v) {
        if (u == v) continue; // self loop
        if (adjMatrix(u, v) == 0) continue; // zero weight edge
        ASSERT(indegree[v] > 0);
        --indegree[v];
        if (indegree[v] == 0) candidates.push_back(v);
      }
    }

    ASSERT(_cached_topord.size() == k, "QG is not acyclic");

    _cached_inverted_topord.clear();
    _cached_inverted_topord.resize(k);
    for (QNodeID u = 0; u < k; ++u) {
      _cached_inverted_topord[_cached_topord[u]] = u;
    }

    _topord_outdated = false;
  }

  const Hypergraph& _hg;
  const Context& _context;
  CycleDetector _cycle_detector;

  std::vector<HypernodeID> _hn_connectivity{};
  std::vector<HypernodeID> _adj_matrix{};
  std::vector<bool> _blocked_moves{};

  bool _topord_outdated{false};
  std::vector<QNodeID> _cached_topord{};
  std::vector<QNodeID> _cached_inverted_topord{};
};
}