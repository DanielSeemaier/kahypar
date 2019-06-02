#include <iostream>
#include <fstream>
#include <sstream>
#include <kahypar/macros.h>

constexpr std::size_t GRAPH_FILENAME = 1;
constexpr std::size_t HGR_FILENAME = 2;
constexpr char COMMENT_CHAR = '#';

int main(int argc, char* argv[]) {
  if (argc != 3) {
    std::cout << "Usage: " << argv[0] << " <graph file> <hgr file>" << std::endl;
    std::exit(-1);
  }

  const std::string graph_filename = argv[GRAPH_FILENAME];
  const std::string hgr_filename = argv[HGR_FILENAME];

  std::ifstream in(graph_filename);
  std::string line;
  bool expect_header = true;

  std::size_t num_nodes = 0;
  std::size_t num_edges = 0;

  const bool has_node_weights = true;
  const bool has_edge_weights = true;

  std::vector<int32_t> hn_weights;
  std::vector<int32_t> he_weights;
  std::vector<std::vector<std::size_t>> hes;

  while (std::getline(in, line)) {
    if (line[0] == COMMENT_CHAR) {
      continue;
    }
    if (expect_header) {
      std::size_t format = 0;
      std::size_t directed = 0;
      std::stringstream(line) >> num_nodes
                              >> num_edges
                              >> format
                              >> directed;
      ASSERT(num_nodes > 0, "Popp graph must contain edges");
      ASSERT(num_edges > 0, "Popp graph must contain edges");
      ASSERT(num_edges <= num_nodes * num_nodes, "Popp graph contains more than n^2 edges");
      ASSERT(format == 11, "Popp graphs must have node and wedge weights");
      ASSERT(directed == 1, "Popp graphs must be directed");
      expect_header = false;
    } else {
      std::stringstream ss(line);

      std::size_t node_weight = 1;
      if (has_node_weights) {
        ss >> node_weight;
      }
      hn_weights.push_back(node_weight);

      std::vector<std::pair<std::size_t, int32_t>> outgoing_edges;
      while (!ss.eof()) {
        std::size_t target;
        std::size_t edge_weight = 1;
        ss >> target;
        target -= 1;
        if (has_edge_weights) {
          ss >> edge_weight;
        }
        outgoing_edges.emplace_back(target, edge_weight);
      }

      ASSERT([&](){
        for (std::size_t i = 0; i + 1 < outgoing_edges.size(); ++i) {
          if (outgoing_edges[i].second != outgoing_edges[i + 1].second) {
            return false;
          }
        }
        return true;
      } (), "Node" << hn_weights.size() << "has outgoing edges with different edge weights");

      if (!outgoing_edges.empty()) {
        he_weights.push_back(outgoing_edges[0].second);

        std::vector<std::size_t> he;
        ASSERT(!hn_weights.empty());
        he.push_back(hn_weights.size() - 1);
        for (const auto& outgoing_edge : outgoing_edges) {
          he.push_back(outgoing_edge.first);
        }
        hes.push_back(std::move(he));
      }
    }
  }
  in.close();

  std::ofstream out(hgr_filename);
  out << "% Generated from " << graph_filename << "\n";
  out << hes.size() << " " << hn_weights.size() << " 11\n";
  for (std::size_t i = 0; i < hes.size(); ++i) {
    ASSERT(he_weights.size() > i);
    out << he_weights[i] << " ";
    for (const auto& pin : hes[i]) {
      out << pin + 1 << " ";
    }
    out << "\n";
  }
  for (const auto& hn_weight : hn_weights) {
    out << hn_weight << "\n";
  }
  out.close();

  return 0;
}