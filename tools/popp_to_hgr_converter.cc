#include <iostream>
#include <fstream>
#include <sstream>

#include "kahypar/macros.h"
#include "kahypar/utils/string.h"

#include "hypergraph_checker.h"

constexpr std::size_t GRAPH_FILENAME = 1;
constexpr std::size_t HGR_FILENAME = 2;
constexpr char COMMENT_CHAR = '#';

enum Output {
  HGR,
  GRAPH
};

using namespace kahypar;

int main(int argc, char* argv[]) {
  if (argc != 3) {
    std::cout << "Usage: " << argv[0] << " <graph file> <*.hgr|*.graph>" << std::endl;
    std::exit(-1);
  }

  const std::string graph_filename = argv[GRAPH_FILENAME];
  const std::string out_filename = argv[HGR_FILENAME];
  Output out_format = (string::ends_with(out_filename, ".graph"))
    ? Output::GRAPH
    : Output::HGR;

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

  if (out_format == Output::HGR) {
    std::ofstream out(out_filename);
    out << "% Generated from " << graph_filename << "\n";
    out << hes.size() << " " << hn_weights.size() << " 11 1\n";
    for (std::size_t i = 0; i < hes.size(); ++i) {
      ASSERT(he_weights.size() > i);
      out << he_weights[i] << " 1 ";
      for (const auto& pin : hes[i]) {
        out << pin + 1 << " ";
      }
      out << "\n";
    }
    for (const auto& hn_weight : hn_weights) {
      out << hn_weight << "\n";
    }
    out.close();
  } else if (out_format == Output::GRAPH) {
    std::ofstream out(out_filename);
    const std::size_t num_edges = std::accumulate(hes.begin(), hes.end(), 0, [](auto acc, const auto& he) {
      ASSERT(he.size() > 1);
      return acc + (he.size() - 1);
    });
    out << hn_weights.size() << " " << num_edges << " 11 1\n";

    ASSERT(hes.size() <= hn_weights.size());
    for (std::size_t hn = 0, he = 0; hn < hn_weights.size(); ++hn) {
      out << hn_weights[hn] << " ";
      if (hes[he][0] == hn) {
        for (std::size_t i = 1; i < hes[he].size(); ++i) {
          out << hes[he][i] + 1 << " " << he_weights[he] << " ";
        }
        ++he;
      }
      out << "\n";
    }
    out.close();
  } else {
    LOG << "Unsupported output format";
    std::exit(-1);
  }

  if (out_format == Output::HGR) {
    validateHypergraphFileAndPanic(out_filename);
  }
  return 0;
}