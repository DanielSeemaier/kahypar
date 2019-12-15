#include "kahypar/io/hypergraph_io.h"

static inline void readBinaryKaffpaD(const std::string& filename,
                                     kaffpa::KaffpaHeader& out_header,
                                     std::vector<kaffpa::Node>& out_nodes,
                                     std::vector<kaffpa::Edge>& out_forward_edges,
                                     std::vector<kaffpa::InEdge>& out_backward_edges,
                                     kaffpa::PartitionConfig& out_config,
                                     std::vector<kaffpa::PartitionID>& out_table,
                                     kaffpa::KaffpaResult& out_result) {
  std::FILE *shm = std::fopen(filename.c_str(), "r");
  if (shm == nullptr) {
    throw std::runtime_error("cannot open binary kaffpaD graph file");
  }

  std::fread(&out_header, sizeof (kaffpa::KaffpaHeader), 1, shm);
  out_nodes.resize(out_header.numberOfNodes + 1);
  out_forward_edges.resize(out_header.numberOfEdges);
  out_backward_edges.resize(out_header.numberOfEdges);
  out_table.resize(out_header.numberOfNodes);

  const auto num_nodes_read = std::fread(out_nodes.data(), sizeof (kaffpa::Node), out_header.numberOfNodes + 1, shm);
  const auto num_forward_edges_read = std::fread(out_forward_edges.data(), sizeof (kaffpa::Edge), out_header.numberOfEdges, shm);
  const auto num_backward_edges_read = std::fread(out_backward_edges.data(), sizeof (kaffpa::InEdge), out_header.numberOfEdges, shm);
  const auto num_configs_read = std::fread(&out_config, sizeof (kaffpa::PartitionConfig), 1, shm);
  const auto num_table_entries_read = std::fread(out_table.data(), sizeof (kaffpa::PartitionID), out_header.numberOfNodes, shm);
  const auto num_results_read = std::fread(&out_result, sizeof (kaffpa::KaffpaResult), 1, shm);

  std::fclose(shm);

  if (num_nodes_read != out_header.numberOfNodes + 1 ||
      num_forward_edges_read != out_header.numberOfEdges ||
      num_backward_edges_read != out_header.numberOfEdges ||
      num_configs_read != 1 ||
      num_table_entries_read != out_header.numberOfNodes ||
      num_results_read != 1) {
    throw std::runtime_error("IO error");
  }
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <*.graph>\n";
    std::exit(1);
  }

  const std::string filename = argv[1];

  kaffpa::KaffpaHeader header{};
  std::vector<kaffpa::Node> nodes;
  std::vector<kaffpa::Edge> forward_edges;
  std::vector<kaffpa::InEdge> backward_edges;
  kaffpa::PartitionConfig config{};
  std::vector<kaffpa::PartitionID> table;
  kaffpa::KaffpaResult result{};

  readBinaryKaffpaD(filename, header, nodes, forward_edges, backward_edges, config, table, result);

  std::cout << "KaffpaResult.k = " << result.k << std::endl;
  std::cout << "KaffpaResult.edgeCut = " << result.edgeCut << std::endl;
  std::cout << "KaffpaResult.objective = " << result.objective << std::endl;

  // compute cut using table
  kaffpa::EdgeWeight cut = 0;
  for (kaffpa::NodeID u = 0; u < header.numberOfNodes; ++u) {
    const auto& node = nodes[u];
    const auto& next = nodes[u + 1];

    const auto part_u = table[u];
    if (part_u >= result.k) {
      std::cout << "Error: invalid block for node " << u << ": " << part_u << std::endl;
      std::cout << "Wont attempt to evaluate the partition table" << std::endl;
      std::exit(1);
    }

    for (kaffpa::EdgeID e = node.firstOutEdge; e < next.firstOutEdge; ++e) {
      const kaffpa::NodeID v = forward_edges[e].target;
      const auto part_v = table[v];
      if (part_u != part_v) {
        cut += forward_edges[e].weight;
      }
    }
  }

  std::cout << "EdgeCut based on the partition table: " << cut << std::endl;
  return 0;
}