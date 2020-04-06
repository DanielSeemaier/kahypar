/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2020 Sebastian Schlag <sebastian.schlag@kit.edu>
 * Copyright (C) 2020 Daniel Seemaier <daniel@seemaier.de>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

namespace kahypar::dhgp {
std::vector<HypernodeID> calculateTopologicalOrdering(const Hypergraph& hg, const bool randomize = true) {
  std::vector<HypernodeID> rank(hg.initialNumNodes());
  for (const HyperedgeID he : hg.edges()) {
    for (const HypernodeID& hh : hg.heads(he)) {
      rank[hh] += hg.edgeNumTails(he);
    }
  }

  std::vector<HypernodeID> candidates;
  for (const HypernodeID hn : hg.nodes()) {
    if (rank[hn] == 0) {
      candidates.push_back(hn);
    }
  }

  std::vector<HypernodeID> topological_ordering;
  topological_ordering.reserve(hg.currentNumNodes());

  while (!candidates.empty()) {
    if (randomize) {
      const auto size = static_cast<int>(candidates.size());
      const auto random_index = Randomize::instance().getRandomInt(0, size - 1);
      std::swap(candidates[random_index], candidates.back());
    }
    const HypernodeID u = candidates.back();
    candidates.pop_back();
    topological_ordering.push_back(u);

    for (const HyperedgeID &he : hg.incidentTailEdges(u)) {
      for (const HypernodeID& hh : hg.heads(he)) {
        ASSERT(rank[hh] > 0);
        --rank[hh];
        if (rank[hh] == 0) candidates.push_back(hh);
      }
    }
  }

  return topological_ordering;
}

std::vector<std::size_t> calculateInvertedTopologicalOrdering(const Hypergraph& hg, const bool randomize = true) {
  const auto topord = calculateTopologicalOrdering(hg, randomize);
  std::vector<std::size_t> inverted_topord(hg.initialNumNodes());
  for (std::size_t i = 0; i < topord.size(); ++i) {
    inverted_topord[topord[i]] = i;
  }
  return inverted_topord;
}

bool checkHypergraphAcyclicity(const Hypergraph& hg) {
  return calculateTopologicalOrdering(hg, false).size() == hg.currentNumNodes();
}

std::vector<HypernodeID> calculateToplevels(const Hypergraph& hg) {
  auto topord = calculateTopologicalOrdering(hg, false);
  std::vector<HypernodeID> toplevels(hg.initialNumNodes());

  for (const HypernodeID& u : topord) {
    for (const HyperedgeID& he : hg.incidentTailEdges(u)) {
      for (const HypernodeID& v : hg.heads(he)) {
        toplevels[v] = std::max(toplevels[v], toplevels[u] + 1);
      }
    }
  }

  return toplevels;
}

std::vector<HypernodeID> calculateReverseToplevels(const Hypergraph& hg) {
  auto topord = calculateTopologicalOrdering(hg, false);
  std::reverse(topord.begin(), topord.end());
  std::vector<HypernodeID> rtoplevels(hg.initialNumNodes());

  HypernodeID max_level = 0;
  for (const HypernodeID& u : topord) {
    for (const HyperedgeID& he : hg.incidentHeadEdges(u)) {
      for (const HypernodeID& v : hg.tails(he)) {
        rtoplevels[v] = std::max(rtoplevels[v], rtoplevels[u] + 1);
        max_level = std::max(max_level, rtoplevels[v]);
      }
    }
  }

  for (auto& rtoplevel : rtoplevels) rtoplevel = max_level - rtoplevel;
  return rtoplevels;
}
} // namespace kahypar::dhgp