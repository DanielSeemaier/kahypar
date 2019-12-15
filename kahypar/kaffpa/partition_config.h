/******************************************************************************
 * partition_config.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#pragma once

#include "kahypar/kaffpa/definitions.h"

namespace kaffpa {

struct KaffpaHeader {
  NodeID numberOfNodes;
  EdgeID numberOfEdges;
  bool hasSecondNodeWeight;
  bool hasNodeFitness;
};

struct KaffpaResult {
  PartitionID k;
  EdgeWeight edgeCut;
  EdgeWeight objective;
  int errorCount;
  int warningCount;
  int numRounds;
};

struct PartitionConfig {
  /* main options */
  unsigned int seed;
  bool kaffpaE;
  bool mpi;

  /* user input */
  NodeWeight min_weight;
  NodeWeight max_weight;
  NodeWeight min_weight2;
  NodeWeight max_weight2;

  /* state */
  bool graph_already_partitioned;
  bool no_new_initial_partitioning;
  bool mh_combine;

  /* set by balance configuration */
  NodeWeight max_node_weight;
  NodeWeight max_node_weight2;
  NodeWeight total_graph_weight;

  /* global search */
  int time_limit;     /* ms */
  int global_partitioning_iterations;
  int global_cycle_iterations;

  /* coarsening */
  stop_rule_t stop_rule;
  int num_vert_stop_factor;

  /* matching */
  matching_type_t matching_type;
  bool edge_rating_tiebreaking;
  edge_rating_t edge_rating;
  permutation_quality_t permutation_quality;
  bool first_level_random_matching;
  bool rate_first_level_inner_outer;
  int aggressive_random_levels;
  bool disable_max_vertex_weight_constraint;
  NodeWeight max_vertex_weight;
  bool gpa_grow_paths_between_blocks;

  /* when using size constraint labeling for clustering */
  node_ordering_type_t node_ordering;
  int cluster_coarsening_factor;
  bool ensemble_clusterings;
  int label_iterations;
  int number_of_clusterings;

  /* initial partitioning */
  bool initial_partitioning_refinement;
  int initial_partitioning_repetitions;
  int initial_partitioning_fill_ratio;

  /* refinement */
  part_refinement_t digraph_refinement;
  bool no_change_convergence;
  permutation_quality_t permutation_during_refinement;
  bool accumulate_on_equal_gain;

  /* evo */
  int mh_local_partitioning_repetitions;
  bool mh_no_mh;
  int mh_flip_coin_mutate_or_combine;
  int mh_flip_coin_mutate_self_or_other;
  int mh_initial_population_fraction;
  bool mh_disable_cross_combine;
  bool mh_cross_combine_original_k;
  bool mh_disable_combine;
  bool mh_enable_quickstart;
  bool mh_diversify;
  bool mh_diversify_best;
  bool mh_enable_tournament_selection;
  bool mh_optimize_communication_volume;
  int mh_pool_size;
  bool mh_optimize_chip;
};

class DefaultPartitionConfig {
public:
  static void set(PartitionConfig &partition_config);
};

} // namespace kaffpa
