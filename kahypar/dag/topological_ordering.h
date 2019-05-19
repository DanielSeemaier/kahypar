#pragma once

#include <vector>
#include <exception>
#include <kahypar/utils/randomize.h>

#include "kahypar/definitions.h"
#include "kahypar/datastructure/hypergraph.h"
#include "kahypar/dag/cycle_detector.h"

namespace kahypar {
namespace dag {
struct CyclicGraphException : public std::exception {
  const char *what() const noexcept {
    return "Hypergraph contains a cycle!";
  }
};

std::vector<HypernodeID> calculateTopologicalOrdering(const Hypergraph& hg) {
  std::vector<HypernodeID> rank(hg.currentNumNodes());
  for (const HyperedgeID& he : hg.edges()) {
    for (const HypernodeID& hh : hg.heads(he)) {
      rank[hh] += hg.edgeNumTails(he);
    }
  }

  std::vector<HypernodeID> candidates;
  for (const HypernodeID& hn : hg.nodes()) {
    if (rank[hn] == 0) {
      candidates.push_back(hn);
    }
  }

  std::vector<HypernodeID> topological_ordering;
  while (!candidates.empty()) {
    HypernodeID u = Randomize::instance().popRandomElement(candidates);
    topological_ordering.push_back(u);

    for (const HyperedgeID &he : hg.incidentEdges(u)) {
      if (hg.isTail(u, he)) {
        for (const HypernodeID& hh : hg.heads(he)) {
          ASSERT(rank[hh] > 0);
          --rank[hh];
          if (rank[hh] == 0) {
            candidates.push_back(hh);
          }
        }
      }
    }
  }

  if (topological_ordering.size() != hg.currentNumNodes()) {
    throw CyclicGraphException{};
  }
  return topological_ordering;
}

bool isAcyclic(const Hypergraph& hg) {
  try {
    calculateTopologicalOrdering(hg);
    return true;
  } catch (const CyclicGraphException& e) {
    return false;
  }
}
} // namespace dag
} // namespace kahypar