/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2014-2016 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include <sys/mman.h>
#include <fcntl.h>

#include "kahypar/application/command_line_options.h"
#include "kahypar/definitions.h"
#include "kahypar/io/hypergraph_io.h"
#include "kahypar/partitioner_facade.h"

#include "kahip_definitions.h"
#include "kahip_config.h"
#include "kahip_misc.h"

using namespace kahypar;

inline Context context(PartitionID k, double epsilon) {
  Context ctx;
  ctx.partition.time_limit = 0;
  ctx.partition.current_v_cycle = 0;
  ctx.partition.graph_partition_filename = "";
  ctx.partition.quiet_mode = true;
  ctx.partition.verbose_output = false;
  ctx.partition.k = k;
  ctx.partition.epsilon = epsilon;
  ctx.partition.final_epsilon = epsilon;
  ctx.partition.mode = Mode::direct_kway;
  ctx.partition.objective = Objective::km1;
  ctx.partition.seed = Randomize::instance().newRandomSeed();
  ctx.partition.hyperedge_size_threshold = 1000;
  ctx.partition.global_search_iterations = 0;
  ctx.preprocessing.enable_min_hash_sparsifier = false;
  ctx.preprocessing.min_hash_sparsifier.min_median_he_size = 28;
  ctx.preprocessing.min_hash_sparsifier.max_hyperedge_size = 1200;
  ctx.preprocessing.min_hash_sparsifier.max_cluster_size = 10;
  ctx.preprocessing.min_hash_sparsifier.min_cluster_size = 2;
  ctx.preprocessing.min_hash_sparsifier.num_hash_functions = 5;
  ctx.preprocessing.min_hash_sparsifier.combined_num_hash_functions = 100;
  ctx.preprocessing.enable_community_detection = false;
  ctx.preprocessing.community_detection.enable_in_initial_partitioning = true;
  ctx.preprocessing.community_detection.reuse_communities = false;
  ctx.preprocessing.community_detection.max_pass_iterations = 100;
  ctx.preprocessing.community_detection.min_eps_improvement = 0.0001;
  ctx.preprocessing.community_detection.edge_weight = LouvainEdgeWeight::hybrid;
  ctx.coarsening.algorithm = CoarseningAlgorithm::ml_style;
  ctx.coarsening.max_allowed_weight_multiplier = 1.0;
  ctx.coarsening.contraction_limit_multiplier = 160;
  ctx.coarsening.rating.rating_function = RatingFunction::heavy_edge;
  ctx.coarsening.rating.community_policy = CommunityPolicy::use_communities;
  ctx.coarsening.rating.heavy_node_penalty_policy = HeavyNodePenaltyPolicy::no_penalty;
  ctx.coarsening.rating.acceptance_policy = AcceptancePolicy::best_prefer_unmatched;
  ctx.coarsening.rating.fixed_vertex_acceptance_policy = FixVertexContractionAcceptancePolicy::fixed_vertex_allowed;
  ctx.initial_partitioning.mode = Mode::recursive_bisection;
  ctx.initial_partitioning.technique = InitialPartitioningTechnique::multilevel;
  ctx.initial_partitioning.coarsening.algorithm = CoarseningAlgorithm::ml_style;
  ctx.initial_partitioning.coarsening.max_allowed_weight_multiplier = 1.0;
  ctx.initial_partitioning.coarsening.contraction_limit_multiplier = 150;
  ctx.initial_partitioning.coarsening.rating.rating_function = RatingFunction::heavy_edge;
  ctx.initial_partitioning.coarsening.rating.community_policy = CommunityPolicy::use_communities;
  ctx.initial_partitioning.coarsening.rating.heavy_node_penalty_policy = HeavyNodePenaltyPolicy::no_penalty;
  ctx.initial_partitioning.coarsening.rating.acceptance_policy = AcceptancePolicy::best_prefer_unmatched;
  ctx.initial_partitioning.coarsening.rating.fixed_vertex_acceptance_policy = FixVertexContractionAcceptancePolicy::fixed_vertex_allowed;
  ctx.initial_partitioning.algo = InitialPartitionerAlgorithm::pool;
  ctx.initial_partitioning.nruns = 20;
  ctx.initial_partitioning.local_search.algorithm = RefinementAlgorithm::twoway_fm;
  ctx.initial_partitioning.local_search.iterations_per_level = -1;
  ctx.initial_partitioning.local_search.fm.stopping_rule = RefinementStoppingRule::simple;
  ctx.initial_partitioning.local_search.fm.max_number_of_fruitless_moves = 50;
  ctx.local_search.algorithm = RefinementAlgorithm::kway_fm_flow_km1;
  ctx.local_search.iterations_per_level = -1;
  ctx.local_search.fm.stopping_rule = RefinementStoppingRule::adaptive_opt;
  ctx.local_search.fm.adaptive_stopping_alpha = 1.0;
  ctx.local_search.fm.max_number_of_fruitless_moves = 350;
  ctx.local_search.flow.algorithm = FlowAlgorithm::ibfs;
  ctx.local_search.flow.alpha = 16;
  ctx.local_search.flow.beta = 128;
  ctx.local_search.flow.network = FlowNetworkType::hybrid;
  ctx.local_search.flow.execution_policy = FlowExecutionMode::exponential;
  ctx.local_search.flow.use_most_balanced_minimum_cut = true;
  ctx.local_search.flow.use_adaptive_alpha_stopping_rule = true;
  ctx.local_search.flow.ignore_small_hyperedge_cut = true;
  ctx.local_search.flow.use_improvement_history = true;
  ctx.partition.max_part_weights.clear();
  return ctx;
}

void readFromShm(const std::string& shm,
                 kahip::KaffpaHeader& header,
                 std::vector<kahip::Node>& nodes,
                 std::vector<kahip::Edge>& forwardEdges,
                 std::vector<kahip::InEdge>& backwardEdges,
                 kahip::PartitionConfig& config) {
  int fd = shm_open(shm.data(), O_RDONLY);
  if (fd < 0) {
    LOG << "Failed to open SHM" << shm << "with RDONLY";
    std::exit(1);
  }

  read(fd, (void*) &header, sizeof(header));

  nodes.reserve(header.numberOfNodes + 1);
  forwardEdges.reserve(header.numberOfEdges);
  backwardEdges.reserve(header.numberOfEdges);

  read(fd, (void*) nodes.data(), sizeof(kahip::Node) * (header.numberOfNodes + 1));
  read(fd, (void*) forwardEdges.data(), sizeof(kahip::Edge) * (header.numberOfEdges));
  read(fd, (void*) backwardEdges.data(), sizeof(kahip::Edge) * (header.numberOfEdges));
  read(fd, (void*) &config, sizeof(kahip::PartitionConfig));

  /*
    mSizeHeader  = sizeof(KaffpaHeader);
    mSizeNodes   = sizeof(Node)        * (mG->getNumNodes() + 1);
    mSizeEdges   = sizeof(Edge)        *  mG->getNumEdges();
    mSizeInEdges = sizeof(InEdge)      *  mG->getNumEdges();
    mSizeConfig  = sizeof(PartitionConfig);
    mSizeTable   = sizeof(PartitionID) *  mG->getNumNodes();
    mSizeResult  = sizeof(KaffpaResult);
   */

  shm_unlink(shm.data());
}

int main(int argc, char* argv[]) {
  if (argc < 3) {
    LOG << "Expected at least 2 arguments:";
    LOG << "\t1st argument: 1=read from shm, 2=read from file";
    LOG << "\t2nd argument: shm name or filename";
    std::exit(-1);
  }

  const int mode = std::stoi(argv[1]);
  const std::string filename = argv[2];

  kahip::KaffpaHeader header;
  std::vector<kahip::Node> nodes;
  std::vector<kahip::Edge> forwardEdges;
  std::vector<kahip::InEdge> backwardEdges;
  kahip::PartitionConfig config;

  if (mode == 1) {
    readFromShm(filename, header, nodes, forwardEdges, backwardEdges, config);
  } else {
    LOG << "TODO";
    std::exit(-1);
  }

  Context ctx = context(config.k, config.epsilon);
  ctx.imbalanced_intermediate_step = false;
  ctx.reduce_balance_during_uncoarsening = false;

  HypernodeID num_hypernodes = header.numberOfNodes;
  HyperedgeID num_hyperedges = header.numberOfEdges;
  HyperedgeIndexVector index_vector;
  HyperedgeVector edge_vector;
  bool is_directed = true;
  std::vector<HypernodeID> num_heads_vector;
  PartitionID num_parts = ctx.partition.k;
  std::vector<HyperedgeWeight> hyperedge_weights;
  std::vector<HypernodeWeight> hypernode_weights;

  index_vector.push_back(edge_vector.size());

  for (HypernodeID hn = 0; hn < header.numberOfNodes; ++hn) {
    for (kahip::EdgeID e = nodes[hn].firstOutEdge; e < nodes[hn + 1].firstOutEdge; ++e) {
      hyperedge_weights.push_back(forwardEdges[e].weight);
      edge_vector.push_back(hn);
      edge_vector.push_back(forwardEdges[e].target);
    }

    index_vector.push_back(edge_vector.size());
    hypernode_weights.push_back(nodes[hn].weight);
    num_heads_vector.push_back(1);
  }

  Hypergraph hypergraph(num_hypernodes, num_hyperedges, index_vector, edge_vector,
                        is_directed, num_heads_vector, num_parts,
                        &hyperedge_weights, &hypernode_weights);

  kahypar::PartitionerFacade().partition(hypergraph, ctx);
  return 0;
}
