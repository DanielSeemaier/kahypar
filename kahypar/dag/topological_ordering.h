#pragma once

#include <vector>
#include <exception>

#include "kahypar/definitions.h"
#include "kahypar/datastructure/hypergraph.h"
#include "kahypar/dag/cycle_detector.h"
#include "kahypar/utils/randomize.h"

namespace kahypar {
namespace dag {
struct CyclicGraphException : public std::exception {
  const char *what() const noexcept override {
    return "Hypergraph contains a cycle!";
  }
};

std::vector<HypernodeID> calculateWeakTopologicalOrdering(const Hypergraph& hg) {
  std::size_t total_deg = 0;

  std::vector<HypernodeID> rank(hg.initialNumNodes());
  for (const HyperedgeID& he : hg.edges()) {
    for (const HypernodeID& hh : hg.heads(he)) {
      rank[hh] += hg.edgeNumTails(he);
      total_deg += hg.edgeNumTails(he);
    }
  }

  std::vector<HypernodeID> todo;
  for (const HypernodeID& hn : hg.nodes()) {
    if (rank[hn] == 0) {
      todo.push_back(hn);
    }
  }

  std::vector<HypernodeID> topological_ordering(hg.initialNumNodes());
  std::size_t level = 0;
  while (!todo.empty()) {
    std::vector<HypernodeID> next_todo;
    for (const HypernodeID& u : todo) {
      topological_ordering[u] = level;

      for (const HyperedgeID &he : hg.incidentEdges(u)) {
        if (hg.isTail(u, he)) {
          for (const HypernodeID& hh : hg.heads(he)) {
            ASSERT(rank[hh] > 0);
            ASSERT(total_deg > 0);
            --rank[hh];
            --total_deg;
            if (rank[hh] == 0) {
              next_todo.push_back(hh);
            }
          }
        }
      }
    }

    ++level;
    todo = next_todo;
  }

  if (total_deg > 0) {
    throw CyclicGraphException{};
  }
  return topological_ordering;
}

std::vector<HypernodeID> calculateTopologicalOrdering(const Hypergraph& hg, bool randomize = true) {
  std::vector<HypernodeID> rank(hg.initialNumNodes());
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
    HypernodeID u;
    if (randomize) {
      u = Randomize::instance().popRandomElement(candidates);
    } else {
      u = candidates.back();
      candidates.pop_back();
    }

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
    calculateTopologicalOrdering(hg, false);
    return true;
  } catch (const CyclicGraphException& e) {
    return false;
  }
}

std::vector<HypernodeID> calculateToplevelValues(const Hypergraph& hg) {
  return calculateWeakTopologicalOrdering(hg);
}
} // namespace dag
} // namespace kahypar