#pragma once

#include <algorithm>
#include <cstdio>

#include "kahypar/partitioner_facade.h"
#include "kahypar/definitions.h"
#include "kahypar/partition/initial_partitioning/i_initial_partitioner.h"
#include "kahypar/partition/initial_partitioning/initial_partitioner_base.h"
#include "kahypar/utils/randomize.h"
#include "kahypar/dag/fix_acyclicity.h"

namespace kahypar {
class KaHyParInitialPartitioner : public IInitialPartitioner, private InitialPartitionerBase<KaHyParInitialPartitioner> {
  using Base = InitialPartitionerBase<KaHyParInitialPartitioner>;
  friend Base;

 public:
  KaHyParInitialPartitioner(Hypergraph& hypergraph, Context& context) :
    Base(hypergraph, context) { }

  ~KaHyParInitialPartitioner() override = default;

  KaHyParInitialPartitioner(const KaHyParInitialPartitioner&) = delete;
  KaHyParInitialPartitioner& operator= (const KaHyParInitialPartitioner&) = delete;

  KaHyParInitialPartitioner(KaHyParInitialPartitioner&&) = delete;
  KaHyParInitialPartitioner& operator= (KaHyParInitialPartitioner&&) = delete;

 private:
  void partitionImpl() final {
    Base::multipleRunsInitialPartitioning();
  }

  void initialPartition() {
    Context ip_context = createContext();
    Hypergraph& hg_copy = _hg;

    const std::string hgr_filename = std::tmpnam(nullptr);
    const std::string part_filename = hgr_filename + ".part" + std::to_string(_context.partition.k) + ".epsilon0.03.seed0.KaHyPar";

    LOG << "Location:" << hgr_filename << part_filename;
    kahypar::io::writeHypergraphFile(hg_copy, hgr_filename, false);

    std::string command = "/Users/danielseemaier/Projects/untouched/kahypar/build/kahypar/application/KaHyPar -h "
      + hgr_filename
      + " -k " + std::to_string(_context.partition.k)
      + " -e 0.03"
      + " -m direct"
      + " -o km1"
      + " -p /Users/danielseemaier/Projects/untouched/kahypar/config/km1_direct_kway_sea17.ini"
      + " --seed 0"
      + " -v true";
    LOG << command;
    std::system(command.c_str());

    std::vector<PartitionID> partition;
    kahypar::io::readPartitionFile(part_filename, partition);
    ASSERT(partition.size() == _hg.initialNumNodes());

//    const auto& map = copy.second;
//    kahypar::PartitionerFacade().partition(hg_copy, ip_context);

//    kahypar::Hypergraph hg_copy(
//      kahypar::io::createHypergraphFromFile("/Users/danielseemaier/Hypergraphs/HGR/ISCAS85/c7552.hgr",
//                                            4));
//    kahypar::PartitionerFacade().partition(hg_copy, ip_context);


    _hg.resetPartitioning();
    _hg.changeK(_context.partition.k);
    for (const HypernodeID& hn : hg_copy.nodes()) {
      //_hg.setNodePart(map[hn], hg_copy.partID(hn));
      _hg.setNodePart(hn, partition[hn]);
    }

    dag::fixAcyclicity(_hg, _context);
    _context.partition.epsilon = metrics::imbalance(_hg, _context);
  }

  Context createContext() const {
    Context ip_context;
    ip_context.partition.current_v_cycle = 0;
    ip_context.partition.graph_partition_filename = "t.part";
    ip_context.partition.quiet_mode = false;
    ip_context.partition.verbose_output = true;
    ip_context.partition.k = _context.partition.k;
    ip_context.partition.epsilon = _context.partition.epsilon;
    ip_context.partition.mode = Mode::direct_kway;
    ip_context.partition.objective = Objective::km1;
    ip_context.partition.seed = -1;
    ip_context.partition.hyperedge_size_threshold = 1000;
    ip_context.partition.global_search_iterations = 0;
    ip_context.preprocessing.enable_min_hash_sparsifier = true;
    ip_context.preprocessing.min_hash_sparsifier.min_median_he_size = 28;
    ip_context.preprocessing.min_hash_sparsifier.max_hyperedge_size = 1200;
    ip_context.preprocessing.min_hash_sparsifier.max_cluster_size = 10;
    ip_context.preprocessing.min_hash_sparsifier.min_cluster_size = 2;
    ip_context.preprocessing.min_hash_sparsifier.num_hash_functions = 5;
    ip_context.preprocessing.min_hash_sparsifier.combined_num_hash_functions = 100;
    ip_context.preprocessing.enable_community_detection = true;
    ip_context.preprocessing.community_detection.enable_in_initial_partitioning = true;
    ip_context.preprocessing.community_detection.reuse_communities = false;
    ip_context.preprocessing.community_detection.max_pass_iterations = 100;
    ip_context.preprocessing.community_detection.min_eps_improvement = 0.0001;
    ip_context.preprocessing.community_detection.edge_weight = LouvainEdgeWeight::hybrid;
    ip_context.coarsening.algorithm = CoarseningAlgorithm::ml_style;
    ip_context.coarsening.max_allowed_weight_multiplier = 1.0;
    ip_context.coarsening.contraction_limit_multiplier = 160;
    ip_context.coarsening.rating.rating_function = RatingFunction::heavy_edge;
    ip_context.coarsening.rating.community_policy = CommunityPolicy::use_communities;
    ip_context.coarsening.rating.heavy_node_penalty_policy = HeavyNodePenaltyPolicy::no_penalty;
    ip_context.coarsening.rating.acceptance_policy = AcceptancePolicy::best_prefer_unmatched;
    ip_context.coarsening.rating.fixed_vertex_acceptance_policy = FixVertexContractionAcceptancePolicy::fixed_vertex_allowed;
    ip_context.initial_partitioning.mode = Mode::recursive_bisection;
    ip_context.initial_partitioning.technique = InitialPartitioningTechnique::multilevel;
    ip_context.initial_partitioning.coarsening.algorithm = CoarseningAlgorithm::ml_style;
    ip_context.initial_partitioning.coarsening.max_allowed_weight_multiplier = 1.0;
    ip_context.initial_partitioning.coarsening.contraction_limit_multiplier = 150;
    ip_context.initial_partitioning.coarsening.rating.rating_function = RatingFunction::heavy_edge;
    ip_context.initial_partitioning.coarsening.rating.community_policy = CommunityPolicy::use_communities;
    ip_context.initial_partitioning.coarsening.rating.heavy_node_penalty_policy = HeavyNodePenaltyPolicy::no_penalty;
    ip_context.initial_partitioning.coarsening.rating.acceptance_policy = AcceptancePolicy::best_prefer_unmatched;
    ip_context.initial_partitioning.coarsening.rating.fixed_vertex_acceptance_policy = FixVertexContractionAcceptancePolicy::fixed_vertex_allowed;
    ip_context.initial_partitioning.algo = InitialPartitionerAlgorithm::pool;
    ip_context.initial_partitioning.nruns = 20;
    ip_context.initial_partitioning.local_search.algorithm = RefinementAlgorithm::twoway_fm;
    ip_context.initial_partitioning.local_search.iterations_per_level = -1;
    ip_context.initial_partitioning.local_search.fm.stopping_rule = RefinementStoppingRule::simple;
    ip_context.initial_partitioning.local_search.fm.max_number_of_fruitless_moves = 50;
    ip_context.local_search.algorithm = RefinementAlgorithm::kway_fm_flow_km1;
    ip_context.local_search.iterations_per_level = -1;
    ip_context.local_search.fm.stopping_rule = RefinementStoppingRule::adaptive_opt;
    ip_context.local_search.fm.adaptive_stopping_alpha = 1.0;
    ip_context.local_search.fm.max_number_of_fruitless_moves = 350;
    ip_context.local_search.flow.algorithm = FlowAlgorithm::ibfs;
    ip_context.local_search.flow.alpha = 16;
    ip_context.local_search.flow.beta = 128;
    ip_context.local_search.flow.network = FlowNetworkType::hybrid;
    ip_context.local_search.flow.execution_policy = FlowExecutionMode::exponential;
    ip_context.local_search.flow.use_most_balanced_minimum_cut = true;
    ip_context.local_search.flow.use_adaptive_alpha_stopping_rule = true;
    ip_context.local_search.flow.ignore_small_hyperedge_cut = true;
    ip_context.local_search.flow.use_improvement_history = true;
    ip_context.partition.max_part_weights.clear();
    return ip_context;
  }
};
} // namespace kahypar