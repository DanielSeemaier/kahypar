#pragma once

#include "kahypar/dag/topological_ordering.h"
#include "kahypar/partition/coarsening/policies/rating_score_policy.h"
#include "kahypar/partition/coarsening/policies/rating_acceptance_policy.h"
#include "kahypar/partition/coarsening/policies/rating_partition_policy.h"

namespace kahypar {
namespace dag {
namespace internal {
bool matchable(const HypernodeID top_u, const HypernodeID top_v) {
  return top_u + 1 == top_v || top_u == top_v + 1;
}
} // namespace internal

std::vector<PartitionID> findAcyclicClustering(const Hypergraph& hg, const Context& context, double max_weight_fraction = 1.0) {
  const auto top = calculateToplevelValues(hg);
  std::vector<PartitionID> leader(hg.initialNumNodes()); // stores the HN that leads the cluster containing the HN
  std::iota(leader.begin(), leader.end(), 0);
  std::vector<HypernodeWeight> weight(hg.initialNumNodes());
  for (const HypernodeID& hn : hg.nodes()) {
    weight[hn] = hg.nodeWeight(hn);
  }
  std::vector<bool> mark(hg.initialNumNodes()); // stores whether a HN is no longer a singleton cluster
  std::vector<HypernodeID> num_bad_neighbors(hg.initialNumNodes());
  std::vector<HypernodeID> leader_bad_neighbors(hg.initialNumNodes(), Hypergraph::kInvalidHypernodeID);

  const auto max_cluster_weight = static_cast<HypernodeWeight>(hg.totalWeight() * max_weight_fraction) + 1;

  std::size_t num_unmatchable = 0;
  std::size_t num_no_neighbor = 0;

  ds::SparseMap<HypernodeID, RatingType> ratings(hg.initialNumNodes(), 0);

  for (const HypernodeID& hn : hg.nodes()) {
    if (mark[hn] || num_bad_neighbors[hn] == 2) {
      if (num_bad_neighbors[hn] == 2) {
        ++num_unmatchable;
      }
      continue;
    }

    HypernodeID best_neighbor = Hypergraph::kInvalidHypernodeID;
    for (const HyperedgeID& he : hg.incidentEdges(hn)) {
      const auto score = HeavyEdgeScore::score(hg, he, context);
      for (const HypernodeID& v : hg.pins(he)) {
        if (v != hn && NormalPartitionPolicy::accept(hg, context, hn, v)) {
          ratings[v] += score;
        }
      }
    }

    // choose best matchable neighbor according to rating
    for (const HyperedgeID& he : hg.incidentEdges(hn)) {
      for (const HypernodeID& neighbor : hg.pins(he)) {
        // hn matchable with every node in neighbor cluster
        bool valid = (!mark[neighbor] && internal::matchable(top[hn], top[leader[neighbor]]))
          || (mark[neighbor] && (top[hn] == top[leader[neighbor]] || top[hn] + 1 == top[leader[neighbor]]));

        // hn matchable with neighbor cluster
        if (valid) {
          valid = num_bad_neighbors[hn] == 0
            || (num_bad_neighbors[hn] == 1 && leader_bad_neighbors[hn] == leader[neighbor]);
        }

        // node pair contractible
        if (valid) {
          valid = hg.partID(hn) == hg.partID(neighbor);
        }

        // cluster doesn't grow too large
        if (valid) {
          valid = weight[leader[neighbor]] + weight[leader[hn]] <= max_cluster_weight;
        }

        if (valid) {
          if (best_neighbor == Hypergraph::kInvalidHypernodeID) {
            best_neighbor = neighbor;
          } else if (ratings[neighbor] > ratings[best_neighbor]
                     || (ratings[neighbor] == ratings[best_neighbor] && Randomize::instance().flipCoin())) {
            best_neighbor = neighbor;
          }
        }
      }
    }
    ratings.clear();

    if (best_neighbor == Hypergraph::kInvalidHypernodeID) {
      ++num_no_neighbor;
      continue;
    }

    HypernodeID l = Hypergraph::kInvalidHypernodeID;
    if (mark[best_neighbor]) {
      l = leader[best_neighbor];
      leader[hn] = l;
      weight[l] += weight[hn];
    } else {
      l = (top[hn] > top[best_neighbor]) ? hn : best_neighbor;
      leader[hn] = l;
      leader[best_neighbor] = l;
      if (l == hn) {
        weight[l] += weight[best_neighbor];
      } else {
        weight[l] += weight[hn];
      }
    }
    ASSERT(l != Hypergraph::kInvalidHypernodeID);

    for (const HyperedgeID& he : hg.incidentEdges(hn)) {
      for (const HypernodeID& neighbor : hg.pins(he)) {
        if (!internal::matchable(top[hn], top[neighbor])) {
          continue;
        }

        if (num_bad_neighbors[neighbor] == 0) {
          num_bad_neighbors[neighbor] = 1;
          leader_bad_neighbors[neighbor] = l;
        } else if (num_bad_neighbors[neighbor] == 1 && leader_bad_neighbors[neighbor] != l) {
          num_bad_neighbors[neighbor] = 2;
        }
      }
    }

    if (!mark[best_neighbor]) {
      for (const HyperedgeID& he : hg.incidentEdges(best_neighbor)) {
        for (const HypernodeID& neighbor : hg.pins(he)) {
          if (!internal::matchable(top[best_neighbor], top[neighbor])) {
            continue;
          }

          if (num_bad_neighbors[neighbor] == 0) {
            num_bad_neighbors[neighbor] = 1;
            leader_bad_neighbors[neighbor] = l;
          } else if (num_bad_neighbors[neighbor] == 1 && leader_bad_neighbors[neighbor] != l) {
            num_bad_neighbors[neighbor] = 2;
          }
        }
      }
      mark[best_neighbor] = true;
    }
    mark[hn] = true;
  }

  ASSERT([&]() {
    std::vector<std::vector<HypernodeID>> clusters(hg.initialNumNodes());
    for (const HypernodeID& hn : hg.nodes()) {
      clusters[leader[hn]].push_back(hn);
    }
    for (const auto& cluster : clusters) {
      if (cluster.empty()) {
        continue;
      }

      HypernodeID level = top[cluster.front()];
      for (const HypernodeID& member : cluster) {
        if (top[member] > level) {
          level = top[member];
        }
      }
      for (const HypernodeID& member : cluster) {
        if (top[member] != level && top[member] + 1 != level) {
          return false;
        }
      }

      HypernodeWeight cluster_weight = 0;
      for (const HypernodeID& member : cluster) {
        cluster_weight += hg.nodeWeight(member);
      }
      if (cluster_weight > max_cluster_weight) {
        return false;
      }
    }

    return true;
  }());

  return leader;
}

void printClusterStats(const Hypergraph& hg, const std::vector<PartitionID>& clustering) {
  const PartitionID max_cluster_num = *std::max_element(clustering.begin(), clustering.end());
  std::vector<std::size_t> cluster_sizes(max_cluster_num + 1);
  std::vector<HypernodeWeight> cluster_weights(max_cluster_num + 1);
  for (const PartitionID& cluster : clustering) {
    ++cluster_sizes[cluster];
  }
  for (const HypernodeID& hn : hg.nodes()) {
    cluster_weights[clustering[hn]] += hg.nodeWeight(hn);
  }

  std::size_t num_clusters = 0;
  std::size_t biggest_cluster = 0;
  std::size_t smallest_cluster = 0;
  HypernodeWeight lightest_cluster = 0;
  HypernodeWeight heaviest_cluster = 0;
  for (std::size_t i = 0; i < cluster_sizes.size(); ++i) {
    if (cluster_sizes[i] > 0) {
      ++num_clusters;
    }
    if (cluster_sizes[i] > cluster_sizes[biggest_cluster]) {
      biggest_cluster = i;
    }
    if (cluster_sizes[i] < cluster_sizes[smallest_cluster] || cluster_sizes[smallest_cluster] == 0) {
      smallest_cluster = i;
    }
    if (cluster_weights[i] < cluster_weights[lightest_cluster] || cluster_weights[lightest_cluster] == 0) {
      lightest_cluster = i;
    }
    if (cluster_weights[i] > cluster_weights[heaviest_cluster]) {
      heaviest_cluster = i;
    }
  }

  LOG << "Number of clusters:" << num_clusters;
  LOG << "Size of biggest cluster:" << cluster_sizes[biggest_cluster];
  LOG << "Size of smallest cluster:" << cluster_sizes[smallest_cluster];
  LOG << "Average cluster size:" << static_cast<double>(clustering.size()) / num_clusters;
  LOG << "Weight of heavierst cluster:" << cluster_weights[heaviest_cluster];
  LOG << "Weight of lightest cluster:" << cluster_weights[lightest_cluster];
  LOG << "Average cluster weight:" << static_cast<double>(hg.totalWeight()) / num_clusters;
}

void printClusterStats(const Hypergraph& hg) {
  ASSERT(!hg.communities().empty(), "No communities!");
  printClusterStats(hg, hg.communities());
}
} // namespace dag
} // namespace kahypar