#pragma once

#include <vector>
#include <exception>
#include <deque>

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
  std::deque<HypernodeID> queue;
  std::vector<HypernodeID> top(hg.initialNumNodes());

  for (const HypernodeID& hn : hg.nodes()) {
    bool source = true;
    for (const HyperedgeID& he : hg.incidentHeadEdges(hn)) {
      for (const HypernodeID& tail : hg.tails(he)) {
        source = false;
        break;
      }
      if (!source) {
        break;
      }
    }
    if (source) {
      queue.push_back(hn);
    }
  }

  while (!queue.empty()) {
    const HypernodeID u = queue.front();
    queue.pop_front();

    for (const HyperedgeID& he : hg.incidentTailEdges(u)) {
      for (const HypernodeID& v : hg.heads(he)) {
        if (top[u] + 1 > top[v]) {
          top[v] = top[u] + 1;
          queue.push_back(v);
        }
      }
    }
  }

  return top;
}
} // namespace dag
} // namespace kahypar