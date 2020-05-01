#pragma once

#ifndef KAHYPAR_ENABLE_DHGP
#error This file should not be included when building without KAHYPAR_ENABLE_DHGP
#endif// KAHYPAR_ENABLE_DHGP

namespace kahypar::dhgp {
template<typename CycleDetector>
class QuotientGraph {
 public:
  QuotientGraph(const Hypergraph& hg, const Context& context)
      : _hg(hg),
        _context(context),
        _cycle_detector(_hg, _context) {}

  void initialize() {
    const PartitionID k = _context.partition.k;
    const HypernodeID n = _hg.initialNumNodes();
    _hn_connectivity.clear();
    _adj_matrix.clear();
    _hn_connectivity.resize(n * k);
    _adj_matrix.resize(k * k);

    for (const HyperedgeID& he : _hg.edges()) {
      for (const HypernodeID& u : _hg.tailPins(he)) {
        for (const HypernodeID& v : _hg.headPins(he)) {
          const PartitionID part_u = _hg.partID(u);
          const PartitionID part_v = _hg.partID(v);

          // u --> v is cut edge
          ++hnConnectivityTo(u, part_v);
          ++hnConnectivityTo(v, part_u);

          if (part_u == part_v) {
            ++adjMatrix(part_u, part_v);
          }
        }
      }
    }

    qgChanged();// flushes blocked moves cache
  }

  bool testAndUpdateBeforeMovement(const HypernodeID hn, const PartitionID to) {
    if (blockedMove(hn, to)) return false;
    return false;
  }

 private:
  HypernodeID& hnConnectivityTo(const HypernodeID hn, const PartitionID part) {
    return _hn_connectivity[hn * _context.partition.k + part];
  }

  HypernodeID& adjMatrix(const HypernodeID u, const HypernodeID v) {
    return _adj_matrix[u * _context.partition.k + v];
  }

  auto& blockedMove(const HypernodeID hn, const PartitionID part) {
    return _blocked_moves[hn * _context.partition.k + part];
  }

  void performMovement(const HypernodeID hn, const PartitionID to) {
    const PartitionID from = _hg.partID(hn);

    for (const HyperedgeID& he : _hg.incidentHeadEdges(hn)) {
      for (const HypernodeID& v : _hg.tailPins(he)) { // hn --> v edge
        const PartitionID v_part = _hg.partID(v);
      }
    }
  }

  void qgChanged() {
    _blocked_moves.clear();
    _blocked_moves.resize(_hg.initialNumNodes() * _context.partition.k);
  }

  const Hypergraph& _hg;
  const Context& _context;
  CycleDetector _cycle_detector;

  std::vector<HypernodeID> _hn_connectivity;
  std::vector<HypernodeID> _adj_matrix;
  std::vector<bool> _blocked_moves;
};
}