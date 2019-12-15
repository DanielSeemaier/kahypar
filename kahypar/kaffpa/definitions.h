/******************************************************************************
 * definitions.h
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

#include <iostream>
#include <limits>
#include <vector>

namespace kaffpa {

/* allows us to disable most of the output during partitioning */
#ifdef DO_CHECK
#define CYCLE_CHECK_OR_FAIL(x, text) if(x) { Log::msg(ERROR, "[cycle_check] " text ": FAIL in " __FILE__ " on %1", __LINE__); Logger::get().printToConsoleAndClear(); exit(1); } \
    else { Log::msg(VERBOSE, "[cycle_check] " text ": OK"); }
#define FULL_CHECK_OR_FAIL(x, text) if(!(x)) { Log::msg(ERROR, "[full_check] " text ": FAIL in " __FILE__ " on %1", __LINE__); Logger::get().printToConsoleAndClear(); exit(1); } \
    else { Log::msg(VERBOSE, "[full_check] " text ": OK"); }
#else
#define CYCLE_CHECK_OR_FAIL(x, text) do {} while (false);
#define FULL_CHECK_OR_FAIL(x, text) do {} while (false);
#endif

typedef int NodeID;
typedef double EdgeRatingType;
typedef int EdgeID;
typedef int InEdgeID;
typedef int PathID;
typedef int PartitionID;
typedef int NodeWeight;
typedef int EdgeWeight;
typedef EdgeWeight Gain;
typedef int Count;

const EdgeID UNDEFINED_EDGE = std::numeric_limits<EdgeID>::max();
const NodeID UNDEFINED_NODE = std::numeric_limits<NodeID>::max();
const PartitionID INVALID_PARTITION = std::numeric_limits<PartitionID>::max();
const int ROOT = 0;

/* matching array has size (no_of_nodes), so for entry in this table we get the matched neighbor */
typedef std::vector<NodeID> coarse_mapping_t;
typedef std::vector<NodeID> matching_t;
typedef std::vector<NodeID> node_permutation_map_t;

/* Coarsening */
typedef enum {
  EXPANSIONSTAR,
  EXPANSIONSTAR2,
  WEIGHT,
  REALWEIGHT,
  PSEUDOGEOM,
  EXPANSIONSTAR2ALGDIST,
  SEPARATOR_MULTX,
  SEPARATOR_ADDX,
  SEPARATOR_MAX,
  SEPARATOR_LOG,
  SEPARATOR_R1,
  SEPARATOR_R2,
  SEPARATOR_R3,
  SEPARATOR_R4,
  SEPARATOR_R5,
  SEPARATOR_R6,
  SEPARATOR_R7,
  SEPARATOR_R8
} edge_rating_t;

typedef enum {
  PERMUTATION_QUALITY_NONE,
  PERMUTATION_QUALITY_FAST,
  PERMUTATION_QUALITY_GOOD
} permutation_quality_t;

typedef enum {
  MATCHING_RANDOM,
  MATCHING_GPA,
  MATCHING_RANDOM_GPA,
  CLUSTER_COARSENING
} matching_type_t;

typedef enum {
  PART_REFINEMENT_NONE = 0,
  PART_REFINEMENT_LOCAL = 1,
  PART_REFINEMENT_GLOBAL = 2,
  PART_REFINEMENT_FM = 3,
  NUM_PART_REFINEMENTS = 4
} part_refinement_t;

typedef enum {
  STOP_RULE_SIMPLE,
  STOP_RULE_MULTIPLE_K,
  STOP_RULE_STRONG
} stop_rule_t;

typedef enum {
  RANDOM_NODEORDERING,
  DEGREE_NODEORDERING
} node_ordering_type_t;

} // namespace kaffpa