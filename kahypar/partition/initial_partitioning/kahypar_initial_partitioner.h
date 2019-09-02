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
    LOG << "initialPartition()";

    _hg.resetPartitioning();
    _hg.changeK(_context.partition.k);
    for (const HypernodeID& hn : _hg.nodes()) {
      _hg.setNodePart(hn, 0);
    }

    performPartition(0, _context.partition.k);

    const bool acyclic = AdjacencyMatrixQuotientGraph<DFSCycleDetector>(_hg, _context).isAcyclic();
    if (!acyclic) {
      LOG << "Error, obtained cyclic IP!";
      std::exit(1);
    } else {
      LOG << "IP is acyclic, nice!";
    }

    for (PartitionID k = 0; k < _context.partition.k; ++k) {
      LOG << "Block" << k << ":" << _hg.partWeight(k) << _hg.partSize(k);
    }

    const double imbalance = metrics::imbalance(_hg, _context);
    if (_context.partition.balance_initial_partition && (imbalance > _context.partition.final_epsilon || imbalance > _context.partition.epsilon)) {
      LOG << "Running hard rebalance to improve IP imbalance from" << imbalance << "to min{"
          << _context.partition.final_epsilon << "," << _context.partition.epsilon << "}";
      rebalancePartition(_hg, _context);
    }
  }

  static void rebalancePartition(Hypergraph& hg, const Context& context, const bool refine_km1 = false) {
    Context balanced_context = context;
    balanced_context.partition.epsilon = context.partition.final_epsilon;
    balanced_context.setupPartWeights(hg.totalWeight());
    LOG << "Improve imbalance to:" << balanced_context.partition.epsilon;

    AdjacencyMatrixQuotientGraph<DFSCycleDetector> qg(hg, balanced_context);
    KMinusOneGainManager gain_manager(hg, balanced_context);
    gain_manager.initialize();
    AcyclicHardRebalanceRefiner hard_balance_refiner(hg, balanced_context, qg, gain_manager);
    hard_balance_refiner.initialize(0);
    Metrics current_metrics = {metrics::hyperedgeCut(hg),
                               metrics::km1(hg),
                               metrics::imbalance(hg, context)};
    UncontractionGainChanges changes{};
    std::vector<HypernodeID> refinement_nodes{};
    hard_balance_refiner.refine(refinement_nodes, {0, 0}, changes, current_metrics);
    //hard_balance_refiner.printSummary();

    if (refine_km1) {
      AcyclicLocalSearchRefiner<NumberOfFruitlessMovesStopsSearch> local_search_refiner(hg, context, qg, gain_manager);
      local_search_refiner.initialize(0);

      std::vector<HypernodeID> local_search_refinement_nodes;
      for (const HypernodeID& hn : hg.nodes()) {
        if (hg.isBorderNode(hn)) {
          local_search_refinement_nodes.push_back(hn);
        }
      }

      local_search_refiner.refine(local_search_refinement_nodes, {0, 0}, changes, current_metrics);
      //local_search_refiner.printSummary();
    }
  }

  void performPartition(const PartitionID part, const PartitionID k) {
    if (k < 2) {
      return;
    }
    LOG << "Subdividing part" << part << "in" << k << "blocks ...";

    const auto pair = extractPartAsUnpartitionedHypergraphForBisection(_hg, part, Objective::km1);
    const auto& hg_ptr = pair.first;
    const auto& map = pair.second;
    hg_ptr->resetPartitioning();
    Context ctx = createContext(2, 0.03);
    kahypar::PartitionerFacade().partition(*hg_ptr, ctx);
    LOG << "Bisection KM1:" << metrics::km1(*hg_ptr);
    LOG << "Bisection imbalance:" << metrics::imbalance(*hg_ptr, ctx);

    PartitionID pre_num_part_0 = 0;
    PartitionID pre_num_part_1 = 0;
    for (const HypernodeID& hn : hg_ptr->nodes()) {
      if (hg_ptr->partID(hn) == 0) {
        ++pre_num_part_0;
      } else {
        ++pre_num_part_1;
      }
    }

    dag::fixBipartitionAcyclicity(*hg_ptr, ctx);
    hg_ptr->initializeNumCutHyperedges();

    LOG << "Bisection KM1 after acyclicity fix:" << metrics::km1(*hg_ptr);
    LOG << "Bisection imbalance after acyclicity fix:" << metrics::imbalance(*hg_ptr, ctx);

    if (_context.initial_partitioning.balance_partition) {
      LOG << "Run HardRebalanceRefiner on initial partition to improve imbalance, then local search to improve KM1";
      rebalancePartition(*hg_ptr, ctx, true);
    }

    LOG << "Bisection KM1 after rebalance + local search:" << metrics::km1(*hg_ptr);
    LOG << "Bisection imbalance after rebalance + local search:" << metrics::imbalance(*hg_ptr, ctx);

    // keep HNs in hg_ptr->partID(hn) == 0 in block `part`
    for (const HypernodeID& hn : hg_ptr->nodes()) {
      if (hg_ptr->partID(hn) == 1) {
        _hg.changeNodePart(map[hn], part, part + 1);
      }
    }

    PartitionID num_part_0 = 0;
    PartitionID num_part_1 = 0;
    for (const HypernodeID& hn : _hg.nodes()) {
      if (_hg.partID(hn) == part) {
        ++num_part_0;
      } else if (_hg.partID(hn) == part + 1) {
        ++num_part_1;
      }
    }

    PartitionID k_part_0 = 0;
    PartitionID k_part_1 = 0;
    if (_context.partition.balance_initial_partition) {
      k_part_0 = k / 2;
      k_part_1 = k / 2;
    } else {
      k_part_0 = std::ceil(k * (static_cast<double>(num_part_0) / hg_ptr->initialNumNodes()));
      k_part_1 = std::floor(k * (static_cast<double>(num_part_1) / hg_ptr->initialNumNodes()));
    }

    k_part_0 = std::max<PartitionID>(1, k_part_0);
    k_part_0 = std::min<PartitionID>(k - 1, k_part_0);
    k_part_1 = std::max<PartitionID>(1, k_part_1);
    k_part_1 = std::min<PartitionID>(k - 1, k_part_1);
    ASSERT(k_part_0 + k_part_1 == k);
    ASSERT(k_part_0 < k); // termination check
    ASSERT(k_part_1 < k); // termination check

    // give up part + 1 in favor of further subdivisions of block `part`
    if (k_part_0 >= 2) {
      for (const HypernodeID& hn : _hg.nodes()) {
        if (_hg.partID(hn) == part + 1) {
          _hg.changeNodePart(hn, part + 1, part + k_part_0);
        }
      }
    } else {
      ASSERT(k_part_0 == 1);
    }

    LOG << "Split" << hg_ptr->initialNumNodes() << "from" << part << "into" << pre_num_part_0 << "and" << pre_num_part_1 << "blocks";
    LOG << "\t\tAfter acyclic fix:" << num_part_0 << "and" << num_part_1 << "blocks";
    LOG << "\tFirst block:" << part << "second block:" << part + k_part_0;
    LOG << "\tContinue with k:" << k_part_0 << "and" << k_part_1;
    const bool acyclic = AdjacencyMatrixQuotientGraph<DFSCycleDetector>(_hg, _context).isAcyclic();
    if (!acyclic) {
      LOG << "Error, obtained cyclic IP!";
      std::exit(1);
    } else {
      LOG << "IP is acyclic, nice!";
    }

    performPartition(part, k_part_0);
    performPartition(part + k_part_0, k_part_1);
  }

  Context createContext(const PartitionID k, const double epsilon) const {
    Context ip_context;
    ip_context.partition.current_v_cycle = 0;
    ip_context.partition.graph_partition_filename = "";
    ip_context.partition.quiet_mode = true;
    ip_context.partition.verbose_output = false;
    ip_context.partition.k = k;
    ip_context.partition.epsilon = epsilon;
    ip_context.partition.final_epsilon = epsilon;
    ip_context.partition.mode = Mode::direct_kway;
    ip_context.partition.objective = Objective::km1;
    ip_context.partition.seed = Randomize::instance().newRandomSeed();
    ip_context.partition.hyperedge_size_threshold = 1000;
    ip_context.partition.global_search_iterations = 0;
    ip_context.preprocessing.enable_min_hash_sparsifier = false;
    ip_context.preprocessing.min_hash_sparsifier.min_median_he_size = 28;
    ip_context.preprocessing.min_hash_sparsifier.max_hyperedge_size = 1200;
    ip_context.preprocessing.min_hash_sparsifier.max_cluster_size = 10;
    ip_context.preprocessing.min_hash_sparsifier.min_cluster_size = 2;
    ip_context.preprocessing.min_hash_sparsifier.num_hash_functions = 5;
    ip_context.preprocessing.min_hash_sparsifier.combined_num_hash_functions = 100;
    ip_context.preprocessing.enable_community_detection = false;
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