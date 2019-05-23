/*
 * Content of this file is mostly copied form hypergraph_io.h
 */
#pragma once

#include <exception>
#include <string>

#include "kahypar/definitions.h"
#include "kahypar/io/hypergraph_io.h"
#include "kahypar/io/partitioning_output.h"
#include "kahypar/partition/metrics.h"
#include "kahypar/dag/topological_ordering.h"
#include "kahypar/utils/math.h"

using namespace std::string_literals;

namespace kahypar {
struct BadHypergraphException : public std::exception {
  explicit BadHypergraphException(std::string message) :
    _message(std::move(message)) {}

  const char *what() const noexcept override {
    return _message.c_str();
  }

 private:
  std::string _message;
};

struct BadHypernodeException : public BadHypergraphException {
  explicit BadHypernodeException(const std::string& message, const HypernodeID hypernode) :
    BadHypergraphException(message + " [hypernode: "s + std::to_string(hypernode) + "]s"),
    _hypernode(hypernode) {}

  HypernodeID hypernode() const { return _hypernode; }

 private:
  HypernodeID _hypernode;
};

struct BadHyperedgeException : public BadHypergraphException {
  explicit BadHyperedgeException(const std::string& message, const HyperedgeID hyperedge) :
    BadHypergraphException(message + " [hyperedge: "s + std::to_string(hyperedge) + "]"s),
    _hyperedge(hyperedge) {}

  HyperedgeID hyperedge() const { return _hyperedge; }

 private:
  HyperedgeID _hyperedge;
};

struct BadPinException : public BadHypergraphException {
  explicit BadPinException(const std::string& message, const HyperedgeID hyperedge, const HypernodeID hypernode) :
    BadHypergraphException(message + " ["s
                           + "hyperedge: "s + std::to_string(hyperedge) + ", "s
                           + "hypernode: "s + std::to_string(hypernode) + "]"s),
    _hyperedge(hyperedge),
    _hypernode(hypernode) {}

  HyperedgeID hyperedge() const { return _hyperedge; }
  HypernodeID hypernode() const { return _hypernode; }

 private:
  HyperedgeID _hyperedge;
  HypernodeID _hypernode;
};

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

void validateHypergraphFile(const std::string& filename) {
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
      throw BadHypergraphException("bad hypergraph type [type: "s + std::to_string(static_cast<int>(hypergraph_type)) + "]"s);
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
        throw BadHyperedgeException("hyperedge is empty", i);
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
          throw BadHypernodeException("invalid hypernode ID", pin + 1);
        }
        if (std::find(edge_vector.begin() + index_vector.back(), edge_vector.end(), pin) != edge_vector.end()) {
          throw BadPinException("hyperedge contains the same pin multiple times", i, pin + 1);
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
    throw BadHypergraphException("hypergraph does not exist [filename: "s + filename + "]"s);
  }
}
} // namespace kahypar