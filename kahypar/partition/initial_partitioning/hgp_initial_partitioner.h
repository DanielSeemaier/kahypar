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
class HgpInitialPartitioner : public IInitialPartitioner, private InitialPartitionerBase<HgpInitialPartitioner> {
  static constexpr bool debug = false;

  using Base = InitialPartitionerBase<HgpInitialPartitioner>;
  friend Base;

 public:
  HgpInitialPartitioner(Hypergraph& hypergraph, Context& context) :
    Base(hypergraph, context) {}

  ~HgpInitialPartitioner() override = default;

  HgpInitialPartitioner(const HgpInitialPartitioner&) = delete;
  HgpInitialPartitioner& operator=(const HgpInitialPartitioner&) = delete;

  HgpInitialPartitioner(HgpInitialPartitioner&&) = delete;
  HgpInitialPartitioner& operator=(HgpInitialPartitioner&&) = delete;

 private:
  void partitionImpl() final {
    Base::multipleRunsInitialPartitioning();
  }

  void initialPartition() {
    DBG << "initialPartition()";
    _total_hgp_time = 0.0;
    HTimer timer;
    timer.start();

    _hg.resetPartitioning();
    _hg.changeK(_context.partition.k);
    for (const HypernodeID& hn : _hg.nodes()) {
      _hg.setNodePart(hn, 0);
    }

    performPartition(0, _context.partition.k);

    const bool acyclic = AdjacencyMatrixQuotientGraph<DFSCycleDetector>(_hg, _context).isAcyclic();
    if (!acyclic) {
      throw std::runtime_error("Produced cyclic initial partition!");
    }

    for (PartitionID k = 0; k < _context.partition.k; ++k) {
      DBG << "Block" << k << ":" << _hg.partWeight(k) << _hg.partSize(k);
    }

    const double imbalance = metrics::imbalance(_hg, _context);
    if (_context.partition.balance_initial_partition &&
        (imbalance > _context.partition.final_epsilon || imbalance > _context.partition.epsilon)) {
      DBG << "Running hard rebalance to improve IP imbalance from" << imbalance << "to min{"
          << _context.partition.final_epsilon << "," << _context.partition.epsilon << "}";
      rebalancePartition(_hg, _context);
    }

    LOG << "Total time spend in blackbox Hypergraph partitioners:" << _total_hgp_time;
    LOG << "Total time spend for initial partitioning:" << timer.stop();
  }

  static void rebalancePartition(Hypergraph& hg, const Context& context, const bool refine_km1 = false) {
    Context balanced_context = context;
    balanced_context.partition.epsilon = context.partition.final_epsilon;
    balanced_context.setupPartWeights(hg.totalWeight());
    DBG << "Improve imbalance to:" << balanced_context.partition.epsilon;

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
      HyperedgeWeight previous_km1 = current_metrics.km1;

      do {
        previous_km1 = current_metrics.km1;
        DBG << "Running KM1 refiner..." << previous_km1;

        AcyclicLocalSearchRefiner<NumberOfFruitlessMovesStopsSearch> local_search_refiner(hg, context, qg,
                                                                                          gain_manager);
        local_search_refiner.initialize(0);
        std::vector<HypernodeID> local_search_refinement_nodes;
        for (const HypernodeID& hn : hg.nodes()) {
          local_search_refinement_nodes.push_back(hn);
        }
        local_search_refiner.refine(local_search_refinement_nodes, {0, 0}, changes, current_metrics);

        DBG << "-->" << current_metrics.km1;
      } while (0.99 * previous_km1 > current_metrics.km1);
    }
  }

//  static void runVCycle(Hypergraph& hg, Context context) {
//    context.setupPartWeights(hg.totalWeight());
//    context.coarsening.algorithm = CoarseningAlgorithm::external;
//    context.coarsening.external_file = "rmlgp";
//    context.coarsening.allow_mixed_contraction = true;
//    context.coarsening.rating.acceptance_policy = AcceptancePolicy::best_prefer_unmatched;
//    context.coarsening.rating.community_policy = CommunityPolicy::ignore_communities;
//    context.coarsening.rating.fixed_vertex_acceptance_policy = FixVertexContractionAcceptancePolicy::fixed_vertex_allowed;
//    context.coarsening.rating.heavy_node_penalty_policy = HeavyNodePenaltyPolicy::no_penalty;
//    context.coarsening.rating.rating_function = RatingFunction::heavy_edge;
//    context.local_search.fm.max_number_of_fruitless_moves = 350;
//    context.local_search.iterations_per_level = 5;
//
//    std::unique_ptr<ICoarsener> coarsener(
//        CoarsenerFactory::getInstance().createObject(
//            context.coarsening.algorithm, hg, context,
//            hg.weightOfHeaviestNode()));
//    AcyclicTwoWayKMinusOneRefiner refiner(hg, context);
//    LOG << context.local_search.iterations_per_level;
//    LOG << "Coarsening with" << context.coarsening.algorithm;
//    LOG << context.local_search.fm.max_number_of_fruitless_moves;
//    coarsener->coarsen(context.coarsening.contraction_limit);
//    LOG << hg.initialNumNodes() << "/" << hg.currentNumNodes();
//    hg.initializeNumCutHyperedges();
//    coarsener->uncoarsen(refiner);
//    refiner.printSummary();
//  }

  void performPartition(const PartitionID part, const PartitionID k) {
    if (k < 2) {
      return;
    }
    DBG << "Subdividing part" << part << "in" << k << "blocks ...";

    const auto pair = extractPartAsUnpartitionedHypergraphForBisection(_hg, part, Objective::km1);
    const auto& hg_ptr = pair.first;
    const auto subgraph_size = hg_ptr->initialNumNodes();
    const auto& map = pair.second;
    hg_ptr->resetPartitioning();
    hg_ptr->setDirected(false);
    invokeHypergraphPartitioner(*hg_ptr, 2, _context.partition.epsilon);
    const bool fix_ip = _context.initial_partitioning.partitioner != "rmlgp";

    Context ctx = createContext(2, _context.partition.epsilon);
    ctx.setupPartWeights(hg_ptr->totalWeight());
    DBG << "Bipartition KM1:" << metrics::km1(*hg_ptr);
    DBG << "Bipartition imbalance:" << metrics::imbalance(*hg_ptr, ctx);

    PartitionID pre_num_part_0 = 0;
    PartitionID pre_num_part_1 = 0;
    for (const HypernodeID& hn : hg_ptr->nodes()) {
      if (hg_ptr->partID(hn) == 0) {
        ++pre_num_part_0;
      } else {
        ++pre_num_part_1;
      }
    }

    ctx.initial_partitioning.balance_partition = _context.initial_partitioning.balance_partition;
    if (fix_ip) {
      dag::fixBipartitionAcyclicity(*hg_ptr, ctx);
    }

//    runVCycle(*hg_ptr, ctx); // TODO refactor

    DBG << "Bipartition KM1 after acyclicity fix + local search:" << metrics::km1(*hg_ptr);
    DBG << "Bipartition imbalance after acyclicity fix + local search:" << metrics::imbalance(*hg_ptr, ctx);
    hg_ptr->printPartSizes();

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
      k_part_0 = std::ceil(k * (static_cast<double>(num_part_0) / subgraph_size));
      k_part_1 = std::floor(k * (static_cast<double>(num_part_1) / subgraph_size));
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

    DBG << "Split" << subgraph_size << "from" << part << "into" << pre_num_part_0 << "and" << pre_num_part_1 << "blocks";
    DBG << "\t\tAfter acyclic fix:" << num_part_0 << "and" << num_part_1 << "blocks";
    DBG << "\tFirst block:" << part << "second block:" << part + k_part_0;
    DBG << "\tContinue with k:" << k_part_0 << "and" << k_part_1;

    performPartition(part, k_part_0);
    performPartition(part + k_part_0, k_part_1);
  }

  void invokeHypergraphPartitioner(Hypergraph& hg, const PartitionID k, const double epsilon) {
    HTimer timer;
    timer.start();
    if (_context.initial_partitioning.partitioner == "kahypar") {
      invokeKaHyPar(hg, k, epsilon);
    } else if (_context.initial_partitioning.partitioner == "patoh") {
      invokePaToH(hg, k, epsilon);
    } else if (_context.initial_partitioning.partitioner == "rmlgp") {
      invokeRMLGP(hg, k, epsilon);
    } else {
      throw std::invalid_argument("bad configuration for i-partitioner");
    }
    _total_hgp_time += timer.stop();
  }

  void invokeRMLGP(Hypergraph& hg, const PartitionID k, const double epsilon) const {
    if (epsilon > 0.03) {
      throw std::runtime_error("unsupported");
    }
    if (k != 2) {
      throw std::runtime_error("unsupported");
    }

    const std::string dot_filename = std::string(std::tmpnam(nullptr)) + ".dot";
    const std::string ip_filename = std::string(std::tmpnam(nullptr));
    io::writeHyperDAGForDotPartitioner(hg, dot_filename);

    std::stringstream cmd_ss;
    cmd_ss << _context.rmlgp_path
      << " " << dot_filename
      << " " << k
      << " --inipart 11"
      << " --ip " << ip_filename
      << " --action ip";
    const std::string cmd = cmd_ss.str();
    LOG << "Calling rMLGP_interface:" << cmd;
    const int rmlgp_exit_code = std::system(cmd.c_str());
    if (rmlgp_exit_code != 0) {
      throw std::runtime_error("rMLGP_interface returned a nonzero exit code");
    }

    std::vector<PartitionID> part;
    io::readPartitionFile(ip_filename, part);
    hg.setPartition(part);

    LOG << "IP with KM1" << metrics::km1(hg) << "imbalance" << metrics::imbalance(hg, _context);
  }

  void invokeKaHyPar(Hypergraph& hg, const PartitionID k, const double epsilon) const {
//    io::writeHypergraphFile(hg, "test.hgr", false);
//    kahypar::Hypergraph hypergraph(
//        kahypar::io::createHypergraphFromFile("test.hgr", _context.partition.k)
//    );
    Context ctx = createContext(k, epsilon);
    kahypar::PartitionerFacade().partition(hg, ctx);
//    hg.resetPartitioning();
//    hg.changeK(k);
//    for (const HypernodeID& hn : hypergraph.nodes()) {
//      hg.setNodePart(hn, hypergraph.partID(hn));
//    }
//    hg.initializeNumCutHyperedges();
  }

  void invokePaToH(Hypergraph& hg, const PartitionID k, const double epsilon) const {
    // write hypergraph to temporary file
    const std::string filename = std::string(std::tmpnam(nullptr));
    io::writeHypergraphForPaToHPartitioning(hg, filename);

    // call PaToH to partition the hypergraph
    std::stringstream command_ss;
    command_ss << _context.patoh_path
      << " " << filename
      << " " << k
      << " UM=O"
      << " SD=" << _context.partition.seed
      << " IB=" << epsilon;
    const std::string command = command_ss.str();
    LOG << "Calling PaToH:" << command;
    const int patoh_exit_code = std::system(command.c_str());
    if (patoh_exit_code != 0) {
      throw std::runtime_error("PaToH returned with a nonzero exit code");
    }

    // make sure PaToH wrote a partition file and load it
    const std::string part_filename = filename + ".part." + std::to_string(k);
    std::ifstream part_file_stream(part_filename);
    if (!part_file_stream.good()) {
      throw std::runtime_error("PaToH didn't create the partition file");
    }
    part_file_stream.close();

    std::vector<PartitionID> part;
    io::readPartitionFile(part_filename, part);
    ASSERT(part.size() == hg.currentNumNodes(), "Mismatch:" << part.size() << "but hypergraph has" << hg.currentNumNodes());

    // apply partition to hypergraph
    hg.resetPartitioning();
    for (const HypernodeID& hn : hg.nodes()) {
      hg.setNodePart(hn, part[hn]);
    }
  }

  static Context createContext(const PartitionID k, const double epsilon) {
    Context ip_context;
    ip_context.partition_evolutionary = false;
    ip_context.partition.time_limit = 0;
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

  double _total_hgp_time{0.0};
};
} // namespace kahypar