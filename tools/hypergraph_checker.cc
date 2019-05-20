/*
 * Content of this file is copied form hypergraph_io.h and slightly modified
 */

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <algorithm>

#include "kahypar/definitions.h"
#include "kahypar/io/hypergraph_io.h"

using namespace kahypar;

static inline void readHGRHeader(std::ifstream& file, HyperedgeID& num_hyperedges,
                                 HypernodeID& num_hypernodes, HypergraphType& hypergraph_type) {
  std::string line;
  std::getline(file, line);

  // skip any comments
  while (line[0] == '%') {
    std::getline(file, line);
  }

  std::istringstream sstream(line);
  int i = 0;
  sstream >> num_hyperedges >> num_hypernodes >> i;
  hypergraph_type = static_cast<HypergraphType>(i);
}

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cout << "Wrong number of arguments!" << std::endl;
    std::cout << "Usage: HypergraphChecker <hypergraph.hgr>" << std::endl;
    return -1;
  }

  std::string filename = argv[1];
  HypernodeID num_hypernodes;
  HyperedgeID num_hyperedges;
  HypergraphType hypergraph_type = HypergraphType::Unweighted;
  HyperedgeIndexVector index_vector;
  HyperedgeVector edge_vector;
  std::ifstream file(filename);

  if (file) {
    readHGRHeader(file, num_hyperedges, num_hypernodes, hypergraph_type);
    if (hypergraph_type != HypergraphType::Unweighted && hypergraph_type != HypergraphType::EdgeWeights
        && hypergraph_type != HypergraphType::NodeWeights && hypergraph_type != HypergraphType::EdgeAndNodeWeights) {
      std::cout << "Bad hypergraph file type" << std::endl;
      std::exit(1);
    }

    const bool has_hyperedge_weights = hypergraph_type == HypergraphType::EdgeWeights ||
                                       hypergraph_type == HypergraphType::EdgeAndNodeWeights;
    const bool has_hypernode_weights = hypergraph_type == HypergraphType::NodeWeights ||
                                       hypergraph_type == HypergraphType::EdgeAndNodeWeights;

    index_vector.reserve(static_cast<size_t>(num_hyperedges) +  /*sentinel*/ 1);
    index_vector.push_back(edge_vector.size());

    std::string line;
    for (HyperedgeID i = 0; i < num_hyperedges; ++i) {
      std::getline(file, line);
      std::istringstream line_stream(line);
      if (line_stream.peek() == EOF) {
        std::cout << "Hyperedge " << i << " is empty" << std::endl;
        std::exit(1);
      }

      if (has_hyperedge_weights) {
        HyperedgeWeight edge_weight;
        line_stream >> edge_weight;
      }
      HypernodeID pin;
      while (line_stream >> pin) {
        // Hypernode IDs start from 0
        --pin;
        if (pin >= num_hypernodes) {
          std::cout << "Invalid hypernode ID: " << pin + 1 << std::endl;
          std::exit(1);
        }
        if (std::find(edge_vector.begin() + index_vector.back(), edge_vector.end(), pin) != edge_vector.end()) {
          std::cout << "Hyperedge " << i << " with multi-pin: " << pin + 1 << std::endl;
          std::exit(1);
        }
        edge_vector.push_back(pin);
      }
      index_vector.push_back(edge_vector.size());
    }

    if (has_hypernode_weights) {
      for (HypernodeID i = 0; i < num_hypernodes; ++i) {
        std::getline(file, line);
        std::istringstream line_stream(line);
        HypernodeWeight node_weight;
        line_stream >> node_weight;
      }
    }
    file.close();
  } else {
    std::cout << "File not found: " << std::endl;
    std::exit(1);
  }

  std::cout << "Hypergraph OK" << std::endl;
  return 0;
}
