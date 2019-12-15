#include <fstream>

#include "kahypar/io/hypergraph_io.h"
#include "hypergraph_checker.h"

enum Format {
  HYPERGRAPH_HGR,
  GRAPH_HGR,
  DGRAPH
};

int main(int argc, char *argv[]) {
  if (argc < 3) {
    std::cout << "Usage: " << argv[0] << " <*.graph> <*.hgr> [hypergraph_hgr | graph_hgr | dgraph]\n";
    std::exit(1);
  }

  const std::string binary_filename = argv[1];
  const std::string out_filename = argv[2];
  Format format = HYPERGRAPH_HGR;
  if (argc > 3) {
    const std::string format_name = argv[3];
    if (format_name == "hypergraph_hgr") {
      format = HYPERGRAPH_HGR;
    } else if (format_name == "graph_hgr") {
      format = GRAPH_HGR;
    } else if (format_name == "dgraph") {
      format = DGRAPH;
    } else {
      throw std::runtime_error("illegal format");
    }
  }

  std::cout << "In: " << binary_filename << "\n";
  std::cout << "Out: " << out_filename << "\n";

  kaffpa::KaffpaHeader header{};
  std::vector<kaffpa::Node> nodes;
  std::vector<kaffpa::Edge> forward_edges;
  std::vector<kaffpa::InEdge> backward_edges;
  kaffpa::PartitionConfig config{};
  kahypar::io::readBinaryKaffpaD(binary_filename, header, nodes, forward_edges, backward_edges, config);
  std::cout << "|V| = " << nodes.size() << ", |E| = " << forward_edges.size() + backward_edges.size() << "\n";

  std::ofstream out(out_filename);

  switch (format) {
    case HYPERGRAPH_HGR: {
      std::cout << "Output format: hypergraph in hypergraph format\n";

      kaffpa::NodeID nonisolated_nodes = 0;
      for (kaffpa::NodeID u = 0; u + 1 < nodes.size(); ++u) {
        const auto &node = nodes[u];
        const auto &next_node = nodes[u + 1];
        if (node.firstOutEdge != next_node.firstOutEdge) {
          ++nonisolated_nodes;
        }
      }

      out << nonisolated_nodes << " " << nodes.size() - 1 <<  " 11 1\n";
      for (kaffpa::NodeID u = 0; u + 1 < nodes.size(); ++u) {
        const auto& node = nodes[u];
        const auto& next_node = nodes[u + 1];
        if (node.firstOutEdge == next_node.firstOutEdge) {
          continue;
        }
        const auto& first_edge = forward_edges[node.firstOutEdge];

        out << first_edge.weight << " 1 " << u + 1 << " ";
        for (kaffpa::EdgeID e = node.firstOutEdge; e < next_node.firstOutEdge; ++e) {
          const auto& edge = forward_edges[e];
          out << edge.target + 1 << " ";
        }
        out << "\n";
      }
      for (kaffpa::NodeID u = 0; u + 1 < nodes.size(); ++u) {
        const auto& node = nodes[u];
        out << node.weight + node.weight2 << "\n";
      }
      break;
    }
    case GRAPH_HGR: {
      std::cout << "Output format: graph in hypergraph format\n";

      out << forward_edges.size() << " " << nodes.size() - 1 << " 11 1\n";
      for (kaffpa::NodeID u = 0; u + 1 < nodes.size(); ++u) {
        const auto& node = nodes[u];
        const auto& next_node = nodes[u + 1];

        for (kaffpa::EdgeID e = node.firstOutEdge; e < next_node.firstOutEdge; ++e) {
          const auto& edge = forward_edges[e];
          out << edge.weight << " 1 " << u + 1 << " " << edge.target + 1 << "\n";
        }
      }
      for (kaffpa::NodeID u = 0; u + 1 < nodes.size(); ++u) {
        const auto& node = nodes[u];
        out << node.weight + node.weight2 << "\n";
      }
      break;
    }
    case DGRAPH: {
      std::cout << "Output format: directed graph in kaffpaD format\n";

      out << nodes.size() - 1 << " " << forward_edges.size() << " 11 1\n";
      for (kaffpa::NodeID u = 0; u + 1 < nodes.size(); ++u) {
        const auto& node = nodes[u];
        const auto& next_node = nodes[u + 1];

        out << node.weight2 + node.weight << " ";
        for (kaffpa::EdgeID e = node.firstOutEdge; e < next_node.firstOutEdge; ++e) {
          const auto& edge = forward_edges[e];
          out << edge.target + 1 << " " << edge.weight << " ";
        }
        out << "\n";
      }
      break;
    }
  }

  out.close();
  if (format == HYPERGRAPH_HGR || format == GRAPH_HGR) {
    std::cout << "Validating ... " << std::flush;
    kahypar::validateHypergraphFile(out_filename);
    std::cout << "OK" << std::endl;
  }

  return 0;
}