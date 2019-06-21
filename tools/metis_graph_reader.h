#pragma once

#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

#include "kahypar/definitions.h"
#include "kahypar/macros.h"

namespace kahypar {
namespace metis {
static constexpr char GRAPH_COMMENT_CHAR = '%';

void readMetisGraph(const std::string& filename,
                    std::vector<int32_t>& node_weights,
                    std::vector<int32_t>& edge_weights,
                    std::vector<std::pair<std::size_t, std::size_t>>& edges) {

  std::ifstream in(filename);
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
}
} // namespace metis
} // namespace kahypar