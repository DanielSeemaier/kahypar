#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

#include "kahypar/partition/metrics.h"
#include "kahypar/partition/context.h"
#include "kahypar/io/hypergraph_io.h"
#include "kahypar/definitions.h"

static constexpr std::size_t HGR_FILENAME = 1;
static constexpr std::size_t GRAPH_FILENAME = 2;
static constexpr std::size_t FIRST_PARTITION_FILENAME = 3;
static constexpr char GRAPH_COMMENT_CHAR = '%';

using namespace kahypar;

static inline double imb(const Hypergraph& hypergraph, const PartitionID k) {
  HypernodeWeight max_weight = hypergraph.partWeight(0);
  for (PartitionID i = 1; i != k; ++i) {
    max_weight = std::max(max_weight, hypergraph.partWeight(i));
  }
  return static_cast<double>(max_weight) /
         ceil(static_cast<double>(hypergraph.totalWeight()) / k) - 1.0;
}

int main(int argc, char* argv[]) {
  if (argc < 4) {
    std::cout << "Usage: " << argv[0] << " <hgr file> <graph file> <partition file>+\n";
    std::exit(-1);
  }

  const std::string hgr_filename = argv[HGR_FILENAME];
  const std::string graph_filename = argv[GRAPH_FILENAME];
  std::vector<std::string> partition_filenames;
  for (std::size_t i = FIRST_PARTITION_FILENAME; i < argc; ++i) {
    partition_filenames.emplace_back(argv[i]);
  }

  // read partition files
  PartitionID max_part = 0;
  std::vector<std::vector<PartitionID>> partitions;
  for (const std::string& partition_filename : partition_filenames) {
    std::vector<PartitionID> partition;
    io::readPartitionFile(partition_filename, partition);
    PartitionID current_max_part = 1;
    for (int block : partition) {
      current_max_part = std::max(current_max_part, block);
    }
    if (max_part != 0 && current_max_part != max_part) {
      LOG << "Error: all partitions must have the same number of blocks";
      std::exit(-1);
    }
    max_part = current_max_part;
    partitions.push_back(std::move(partition));
  }

  // read hypergraph
  // TODO integrate this into the file format
  constexpr bool directed = true;
  constexpr HypernodeID num_heads_per_edge = 1;
  Hypergraph hypergraph(io::createHypergraphFromFile(hgr_filename, max_part + 1, directed, num_heads_per_edge));

  // read graph
  std::vector<int32_t> node_weights;
  std::vector<int32_t> edge_weights;
  std::vector<std::pair<std::size_t, std::size_t>> edges;

  std::ifstream in(graph_filename);
  std::size_t num_nodes = 0;
  std::size_t num_edges = 0;
  bool has_node_weights = false;
  bool has_edge_weights = false;
  std::string line;
  bool expect_header = true;
  while (std::getline(in, line)) {
    if (line[0] == GRAPH_COMMENT_CHAR) {
      continue;
    }

    if (expect_header) {
      std::size_t format;
      std::stringstream(line) >> num_nodes
                              >> num_edges
                              >> format;
      has_node_weights = (format == 10 || format == 11);
      has_edge_weights = (format == 1 || format == 11);
      expect_header = false;
    } else {
      std::stringstream ss(line);
      std::size_t node_weight = 1;
      if (has_node_weights) {
        ss >> node_weight;
      }
      node_weights.push_back(node_weight);
      ASSERT(hypergraph.currentNumNodes() >= node_weights.size());
      if (hypergraph.nodeWeight(node_weights.size() - 1) != node_weight) {
        LOG << "Error: node" << node_weights.size() - 1
            << "has node weight" << node_weight
            << "in graph file but" << hypergraph.nodeWeight(node_weights.size() - 1)
            << "in hypergraph file";
        std::exit(-1);
      }

      std::size_t head = 0;
      while (ss >> head) {
        std::size_t edge_weight = 1;
        if (has_edge_weights) {
          ss >> edge_weight;
        }
        edge_weights.push_back(edge_weight);
        ASSERT(head > 0);
        ASSERT(!node_weights.empty());
        edges.emplace_back(node_weights.size() - 1, head - 1);
      }
    }
  }
  in.close();
  ASSERT(edge_weights.size() == edges.size());

  if (node_weights.size() != num_nodes) {
    LOG << "Error: bad number of nodes in graph file";
    std::exit(-1);
  }
  if (edges.size() != num_edges) {
    LOG << "Error: bad number of edges in graph file";
    std::exit(-1);
  }
  if (num_nodes != hypergraph.currentNumNodes()) {
    LOG << "Error: graph contains" << num_nodes << "nodes but hypergraph contains" << hypergraph.currentNumNodes()
        << "nodes";
    std::exit(-1);
  }

  LOG << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%";
  LOG << "% k: number of blocks                                                          %";
  LOG << "% imbalance: max_weight / total_weight - 1.0                                   %";
  LOG << "% km1: sum over all edges e: w(e) * (lambda - 1) [Hypergraph metric]           %";
  LOG << "% cut: sum over cut edges e: w(e) [Graph metric]                               %";
  LOG << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%";

  for (std::size_t i = 0; i < partition_filenames.size(); ++i) {
    const auto& partition_filename = partition_filenames[i];
    const auto& partition = partitions[i];

    hypergraph.resetPartitioning();
    for (size_t index = 0; index < partition.size(); ++index) {
      hypergraph.setNodePart(index, partition[index]);
    }
    Context context;
    context.partition.k = max_part + 1;

    std::cout << "k(" << partition_filename << ")=" << max_part + 1 << std::endl;
    std::cout << "imbalance(" << partition_filename << ")=" << imb(hypergraph, context.partition.k) << std::endl;
    std::cout << "km1(" << partition_filename << ")=" << metrics::km1(hypergraph) << std::endl;

    int32_t graph_cut = 0;
    for (std::size_t j = 0; j < edges.size(); ++j) {
      const auto& weight = edge_weights[j];
      const auto& tail = edges[j].first;
      const auto& head = edges[j].second;
      if (partition[tail] != partition[head]) {
        graph_cut += weight;
      }
    }

    std::cout << "cut(" << partition_filename << ")=" << graph_cut << std::endl;
    LOG << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%";
  }

  return 0;
}