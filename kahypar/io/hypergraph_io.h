/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2014-2016 Sebastian Schlag <sebastian.schlag@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "kahypar/definitions.h"
#include "kahypar/kaffpa/partition_config.h"
#include "kahypar/kaffpa/definitions.h"
#include "kahypar/kaffpa/graph_access.h"

namespace kahypar {
namespace io {
using Mapping = std::unordered_map<HypernodeID, HypernodeID>;

static inline void readHGRHeader(std::ifstream &file, HyperedgeID &num_hyperedges,
                                 HypernodeID &num_hypernodes, HypergraphType &hypergraph_type,
                                 bool &is_directed) {
  std::string line;
  std::getline(file, line);

  // skip any comments
  while (line[0] == '%') {
    std::getline(file, line);
  }

  std::istringstream sstream(line);
  int i = 0;
  sstream >> num_hyperedges >> num_hypernodes >> i >> is_directed;
  hypergraph_type = static_cast<HypergraphType>(i);
}

static inline void readHypergraphFile(const std::string &filename, HypernodeID &num_hypernodes,
                                      HyperedgeID &num_hyperedges,
                                      HyperedgeIndexVector &index_vector,
                                      HyperedgeVector &edge_vector,
                                      bool &is_directed,
                                      NumHeadsVector &num_heads_vector,
                                      HyperedgeWeightVector *hyperedge_weights = nullptr,
                                      HypernodeWeightVector *hypernode_weights = nullptr) {
  ASSERT(!filename.empty(), "No filename for hypergraph file specified");
  HypergraphType hypergraph_type = HypergraphType::Unweighted;
  std::ifstream file(filename);
  if (file) {
    readHGRHeader(file, num_hyperedges, num_hypernodes, hypergraph_type, is_directed);
    ASSERT(hypergraph_type == HypergraphType::Unweighted ||
           hypergraph_type == HypergraphType::EdgeWeights ||
           hypergraph_type == HypergraphType::NodeWeights ||
           hypergraph_type == HypergraphType::EdgeAndNodeWeights,
           "Hypergraph in file has wrong type");

    const bool has_hyperedge_weights = hypergraph_type == HypergraphType::EdgeWeights ||
                                       hypergraph_type == HypergraphType::EdgeAndNodeWeights ?
                                       true : false;
    const bool has_hypernode_weights = hypergraph_type == HypergraphType::NodeWeights ||
                                       hypergraph_type == HypergraphType::EdgeAndNodeWeights ?
                                       true : false;

    index_vector.reserve(static_cast<size_t>(num_hyperedges) +  /*sentinel*/ 1);
    index_vector.push_back(edge_vector.size());

    std::string line;
    for (HyperedgeID i = 0; i < num_hyperedges; ++i) {
      std::getline(file, line);
      std::istringstream line_stream(line);
      if (line_stream.peek() == EOF) {
        std::cerr << "Error: Hyperedge " << i << " is empty" << std::endl;
        exit(1);
      }

      if (has_hyperedge_weights) {
        HyperedgeWeight edge_weight;
        line_stream >> edge_weight;
        if (hyperedge_weights == nullptr) {
          LOG << "****** ignoring hyperedge weights ******";
        } else {
          ASSERT(hyperedge_weights != nullptr, "Hypergraph has hyperedge weights");
          hyperedge_weights->push_back(edge_weight);
        }
      }
      HypernodeID num_heads = 0;
      if (is_directed) {
        line_stream >> num_heads;
      }
      num_heads_vector.push_back(num_heads);
      HypernodeID pin;
      while (line_stream >> pin) {
        // Hypernode IDs start from 0
        --pin;
        ASSERT(pin < num_hypernodes, "Invalid hypernode ID");
        edge_vector.push_back(pin);
      }
      index_vector.push_back(edge_vector.size());
    }

    if (has_hypernode_weights) {
      if (hypernode_weights == nullptr) {
        LOG << " ****** ignoring hypernode weights ******";
      } else {
        ASSERT(hypernode_weights != nullptr, "Hypergraph has hypernode weights");
        for (HypernodeID i = 0; i < num_hypernodes; ++i) {
          std::getline(file, line);
          std::istringstream line_stream(line);
          HypernodeWeight node_weight;
          line_stream >> node_weight;
          hypernode_weights->push_back(node_weight);
        }
      }
    }
    file.close();
  } else {
    std::cerr << "Error: File not found: " << std::endl;
  }
}

static inline void readHypergraphFile(const std::string &filename,
                                      HypernodeID &num_hypernodes,
                                      HyperedgeID &num_hyperedges,
                                      std::unique_ptr<size_t[]> &index_vector,
                                      std::unique_ptr<HypernodeID[]> &edge_vector,
                                      bool &is_directed,
                                      std::unique_ptr<HypernodeID[]> &num_heads_vector,
                                      std::unique_ptr<HyperedgeWeight[]> &hyperedge_weights,
                                      std::unique_ptr<HypernodeWeight[]> &hypernode_weights) {
  HyperedgeIndexVector index_vec;
  HyperedgeVector edge_vec;
  NumHeadsVector num_heads_vec;
  HyperedgeWeightVector edge_weights_vec;
  HypernodeWeightVector node_weights_vec;

  readHypergraphFile(filename, num_hypernodes, num_hyperedges, index_vec,
                     edge_vec, is_directed, num_heads_vec, &edge_weights_vec, &node_weights_vec);

  ASSERT(index_vector == nullptr);
  ASSERT(edge_vector == nullptr);
  ASSERT(num_heads_vector == nullptr);
  index_vector = std::make_unique<size_t[]>(index_vec.size());
  edge_vector = std::make_unique<HypernodeID[]>(edge_vec.size());
  num_heads_vector = std::make_unique<HypernodeID[]>(num_heads_vec.size());

  memcpy(index_vector.get(), index_vec.data(), index_vec.size() * sizeof(size_t));
  memcpy(edge_vector.get(), edge_vec.data(), edge_vec.size() * sizeof(HypernodeID));
  memcpy(num_heads_vector.get(), num_heads_vec.data(), num_heads_vec.size() * sizeof(HypernodeID));

  if (!edge_weights_vec.empty()) {
    ASSERT(hyperedge_weights == nullptr);
    hyperedge_weights = std::make_unique<HyperedgeWeight[]>(edge_weights_vec.size());
    memcpy(hyperedge_weights.get(), edge_weights_vec.data(),
           edge_weights_vec.size() * sizeof(HyperedgeWeight));
  }

  if (!node_weights_vec.empty()) {
    ASSERT(hypernode_weights == nullptr);
    hypernode_weights = std::make_unique<HypernodeWeight[]>(node_weights_vec.size());
    memcpy(hypernode_weights.get(), node_weights_vec.data(),
           node_weights_vec.size() * sizeof(HypernodeWeight));
  }
}

static inline Hypergraph createHypergraphFromFile(const std::string &filename,
                                                  const PartitionID num_parts) {
  HypernodeID num_hypernodes;
  HyperedgeID num_hyperedges;
  HyperedgeIndexVector index_vector;
  HyperedgeVector edge_vector;
  HypernodeWeightVector hypernode_weights;
  HyperedgeWeightVector hyperedge_weights;
  NumHeadsVector num_heads_vector;
  bool is_directed;
  readHypergraphFile(filename, num_hypernodes, num_hyperedges,
                     index_vector, edge_vector, is_directed, num_heads_vector,
                     &hyperedge_weights, &hypernode_weights);
  return Hypergraph(num_hypernodes, num_hyperedges, index_vector, edge_vector,
                    is_directed, num_heads_vector, num_parts,
                    &hyperedge_weights, &hypernode_weights);
}


static inline void writeHypernodeWeights(std::ofstream &out_stream, const Hypergraph &hypergraph) {
  for (const HypernodeID &hn : hypergraph.nodes()) {
    out_stream << hypergraph.nodeWeight(hn) << "\n";
  }
}

static inline void
writeHGRHeader(std::ofstream &out_stream, const Hypergraph &hypergraph, const bool respect_directed = true) {
  out_stream << hypergraph.initialNumEdges() << " " << hypergraph.initialNumNodes() << " ";
  if (hypergraph.type() != HypergraphType::Unweighted) {
    out_stream << static_cast<int>(hypergraph.type());
  }
  if (respect_directed && hypergraph.isDirected()) {
    out_stream << " 1";
  }
  out_stream << std::endl;
}

static inline void
writeHypergraphFile(const Hypergraph &hypergraph, const std::string &filename, const bool respect_directed = true) {
  ASSERT(!filename.empty(), "No filename for hypergraph file specified");
  ALWAYS_ASSERT(!hypergraph.isModified(), "Hypergraph is modified. Reindexing HNs/HEs necessary.");

  std::ofstream out_stream(filename.c_str());
  writeHGRHeader(out_stream, hypergraph, respect_directed);

  for (const HyperedgeID &he : hypergraph.edges()) {
    if (hypergraph.type() == HypergraphType::EdgeWeights ||
        hypergraph.type() == HypergraphType::EdgeAndNodeWeights) {
      out_stream << hypergraph.edgeWeight(he) << " ";
    }
    if (respect_directed && hypergraph.isDirected()) {
      out_stream << hypergraph.edgeNumHeads(he) << " ";
    }
    for (const HypernodeID &pin : hypergraph.pins(he)) {
      out_stream << pin + 1 << " ";
    }
    out_stream << "\n";
  }

  if (hypergraph.type() == HypergraphType::NodeWeights ||
      hypergraph.type() == HypergraphType::EdgeAndNodeWeights) {
    writeHypernodeWeights(out_stream, hypergraph);
  }
  out_stream.close();
}


static inline void writeHypergraphToGraphMLFile(const Hypergraph &hypergraph,
                                                const std::string &filename,
                                                const std::vector<PartitionID> *hn_cluster_ids = nullptr,
                                                const std::vector<PartitionID> *he_cluster_ids = nullptr) {
  std::ofstream out_stream(filename.c_str());

  out_stream << R"(<?xml version="1.0" encoding="UTF-8" standalone="no"?>)"
             << R"( <graphml xmlns="http://graphml.graphdrawing.org/xmlns")"
             << R"( xmlns:java="http://www.yworks.com/xml/yfiles-common/1.0/java")"
             << R"( xmlns:sys="http://www.yworks.com/xml/yfiles-common/markup/primitives/2.0")"
             << R"( xmlns:x="http://www.yworks.com/xml/yfiles-common/markup/2.0")"
             << R"( xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance")"
             << R"( xmlns:y="http://www.yworks.com/xml/graphml")"
             << R"( xmlns:yed="http://www.yworks.com/xml/yed/3")"
             << R"( xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns)"
             << R"(http://www.yworks.com/xml/schema/graphml/1.1/ygraphml.xsd">)"
             << std::endl;

  out_stream << R"(<key id="d0" for="node" attr.name="weight" attr.type="double"/>)" << std::endl;
  out_stream << R"(<key id="d1" for="node" attr.name="part" attr.type="int"/>)" << std::endl;
  out_stream << R"(<key id="d2" for="node" attr.name="iscutedge" attr.type="int"/>)" << std::endl;
  out_stream << R"(<key id="d7" for="node" attr.name="modclass" attr.type="int"/>)" << std::endl;
  out_stream << R"(<key id="d8" for="node" attr.name="color" attr.type="string"/>)" << std::endl;
  out_stream << R"(<graph id="G" edgedefault="undirected">)" << std::endl;
  for (const HypernodeID &hn : hypergraph.nodes()) {
    out_stream << R"(<node id="n)" << hn << R"(">)" << std::endl;
    out_stream << R"(<data key="d0">)" << hypergraph.nodeWeight(hn) << "</data>" << std::endl;
    if (hn_cluster_ids != nullptr) {
      out_stream << R"(<data key="d7">)" << (*hn_cluster_ids)[hn] << "</data>" << std::endl;
    } else {
      out_stream << R"(<data key="d1">)" << hypergraph.partID(hn) << "</data>" << std::endl;
    }

    out_stream << R"(<data key="d2">)" << 42 << "</data>" << std::endl;
    out_stream << R"(<data key="d8">)" << "blue" << "</data>" << std::endl;
    out_stream << "</node>" << std::endl;
  }

  HyperedgeID edge_id = 0;
  for (const HyperedgeID &he : hypergraph.edges()) {
    // const HyperedgeID he_id = hypergraph.initialNumNodes() + he;
    out_stream << R"(<node id="h)" << he << R"(">)" << std::endl;
    out_stream << R"(<data key="d0">)" << hypergraph.edgeWeight(he) << "</data>" << std::endl;
    if (he_cluster_ids != nullptr) {
      out_stream << R"(<data key="d7">)" << (*he_cluster_ids)[he] << "</data>" << std::endl;
    } else {
      out_stream << R"(<data key="d1">)" << -1 << "</data>" << std::endl;
    }
    out_stream << R"(<data key="d2">)" << (hypergraph.connectivity(he) > 1) << "</data>" << std::endl;
    out_stream << R"(<data key="d8">)" << "red" << "</data>" << std::endl;
    out_stream << "</node>" << std::endl;
    for (const HypernodeID &pin : hypergraph.pins(he)) {
      out_stream << R"(<edge id="e)" << edge_id++ << R"(" source="n)" << pin << R"(" target="h)"
                 << he << R"("/>)" << std::endl;
    }
  }

  out_stream << "</graph>" << std::endl;
  out_stream << "</graphml>" << std::endl;
  out_stream.close();
}


static inline void writeHypergraphForhMetisPartitioning(const Hypergraph &hypergraph,
                                                        const std::string &filename,
                                                        const Mapping &mapping) {
  ASSERT(!filename.empty(), "No filename for hMetis initial partitioning file specified");
  std::ofstream out_stream(filename.c_str());

  // coarse graphs always have edge and node weights, even if graph wasn't coarsend
  out_stream << hypergraph.currentNumEdges() << " " << hypergraph.currentNumNodes() << " ";
  out_stream << static_cast<int>(HypergraphType::EdgeAndNodeWeights);
  out_stream << std::endl;

  for (const HyperedgeID &he : hypergraph.edges()) {
    out_stream << hypergraph.edgeWeight(he) << " ";
    for (const HypernodeID &pin : hypergraph.pins(he)) {
      ASSERT(mapping.find(pin) != mapping.end(), "No mapping found for pin " << pin);
      out_stream << mapping.find(pin)->second + 1 << " ";
    }
    out_stream << std::endl;
  }

  writeHypernodeWeights(out_stream, hypergraph);
  out_stream.close();
}

static inline void writeHypergraphForPaToHPartitioning(const Hypergraph &hypergraph,
                                                       const std::string &filename,
                                                       const Mapping &mapping) {
  ASSERT(!filename.empty(), "No filename for PaToH initial partitioning file specified");
  std::ofstream out_stream(filename.c_str());
  out_stream << 1;                     // 1-based indexing
  out_stream << " " << hypergraph.currentNumNodes() << " " << hypergraph.currentNumEdges() << " "
             << hypergraph.currentNumPins();
  out_stream << " " << 3 << std::endl;  // weighting scheme: both edge and node weights

  for (const HyperedgeID &he : hypergraph.edges()) {
    out_stream << hypergraph.edgeWeight(he) << " ";
    for (const HypernodeID &pin : hypergraph.pins(he)) {
      ASSERT(mapping.find(pin) != mapping.end(), "No mapping found for pin " << pin);
      out_stream << mapping.find(pin)->second + 1 << " ";
    }
    out_stream << std::endl;
  }

  for (const HypernodeID &hn : hypergraph.nodes()) {
    out_stream << hypergraph.nodeWeight(hn) << " ";
  }
  out_stream << std::endl;
  out_stream.close();
}

static inline void writeHypergraphForPaToHPartitioning(const Hypergraph &hypergraph,
                                                       const std::string &filename) {
  ASSERT(!filename.empty(), "No filename for PaToH initial partitioning file specified");
  std::ofstream out_stream(filename.c_str());
  out_stream << 0;                     // 0-based indexing
  out_stream << " " << hypergraph.currentNumNodes() << " " << hypergraph.currentNumEdges() << " "
             << hypergraph.currentNumPins();
  out_stream << " " << 3 << std::endl;  // weighting scheme: both edge and node weights

  for (const HyperedgeID &he : hypergraph.edges()) {
    out_stream << hypergraph.edgeWeight(he) << " ";
    for (const HypernodeID &pin : hypergraph.pins(he)) {
      // ASSERT(mapping.find(pin) != mapping.end(), "No mapping found for pin " << pin);
      out_stream << pin << " ";
    }
    out_stream << std::endl;
  }

  for (const HypernodeID &hn : hypergraph.nodes()) {
    out_stream << hypergraph.nodeWeight(hn) << " ";
  }
  out_stream << std::endl;
  out_stream.close();
}

static inline void writeHyperDAGForDotPartitioner(const Hypergraph &hypergraph, const std::string &filename) {
  std::ofstream dot(filename);
  dot << "digraph G {\n";
  for (const HypernodeID &hn : hypergraph.nodes()) {
    dot << hn << ";\n";
  }
  for (const HyperedgeID &he : hypergraph.edges()) {
    for (const HypernodeID &tail : hypergraph.tails(he)) {
      for (const HypernodeID &head : hypergraph.heads(he)) {
        dot << tail << "->" << head << ";\n";
      }
    }
  }
  dot << "}\n";
}


static inline void readPartitionFile(const std::string &filename, std::vector<PartitionID> &partition) {
  ASSERT(!filename.empty(), "No filename for partition file specified");
  ASSERT(partition.empty(), "Partition vector is not empty");
  std::ifstream file(filename);
  if (file) {
    int part;
    while (file >> part) {
      partition.push_back(part);
    }
    file.close();
  } else {
    std::cerr << "Error: File not found: " << std::endl;
  }
}

static inline void writePartitionFile(const Hypergraph &hypergraph, const std::string &filename) {
  if (filename.empty()) {
    LOG << "No filename for partition file specified";
  } else {
    std::ofstream out_stream(filename.c_str());
    for (const HypernodeID &hn : hypergraph.nodes()) {
      out_stream << hypergraph.partID(hn) << std::endl;
    }
    out_stream.close();
  }
}

static inline void readFixedVertexFile(Hypergraph &hypergraph, const std::string &filename) {
  ASSERT(!filename.empty(), "No filename for partition file specified");
  std::ifstream file(filename);
  if (file) {
    PartitionID part;
    HypernodeID hn = 0;
    while (file >> part) {
      if (part != -1) {
        hypergraph.setFixedVertex(hn, part);
      }
      hn++;
    }
    file.close();
  } else {
    std::cerr << "Error: File not found: " << filename << std::endl;
  }
}

static inline void writeFixedVertexFile(const Hypergraph &hypergraph, const std::string &filename) {
  ASSERT(!filename.empty(), "No filename for partition file specified");
  std::ofstream out_stream(filename.c_str());
  for (const HypernodeID &hn : hypergraph.nodes()) {
    out_stream << hypergraph.fixedVertexPartID(hn) << std::endl;
  }
  out_stream.close();
}

static inline void readBinaryKaffpaD(const std::string &filename,
                                     kaffpa::KaffpaHeader &out_header,
                                     std::vector<kaffpa::Node> &out_nodes,
                                     std::vector<kaffpa::Edge> &out_forward_edges,
                                     std::vector<kaffpa::InEdge> &out_backward_edges,
                                     kaffpa::PartitionConfig &out_config) {
  std::FILE *shm = std::fopen(filename.c_str(), "r");
  if (shm == nullptr) {
    throw std::runtime_error("cannot open binary kaffpaD graph file");
  }

  std::fread(&out_header, sizeof(kaffpa::KaffpaHeader), 1, shm);
  out_nodes.resize(out_header.numberOfNodes + 1);
  out_forward_edges.resize(out_header.numberOfEdges);
  out_backward_edges.resize(out_header.numberOfEdges);

  const auto num_nodes_read = std::fread(out_nodes.data(), sizeof(kaffpa::Node), out_header.numberOfNodes + 1, shm);
  const auto num_forward_edges_read = std::fread(out_forward_edges.data(), sizeof(kaffpa::Edge),
                                                 out_header.numberOfEdges, shm);
  const auto num_backward_edges_read = std::fread(out_backward_edges.data(), sizeof(kaffpa::InEdge),
                                                  out_header.numberOfEdges, shm);
  const auto num_configs_read = std::fread(&out_config, sizeof(kaffpa::PartitionConfig), 1, shm);
  std::fclose(shm);

  if (num_nodes_read != out_header.numberOfNodes + 1 ||
      num_forward_edges_read != out_header.numberOfEdges ||
      num_backward_edges_read != out_header.numberOfEdges ||
      num_configs_read != 1) {
    throw std::runtime_error("IO error");
  }
}

static inline void buildBinaryKaffpaDHNMap(const std::vector<kaffpa::Node> &nodes,
                                           const std::vector<kaffpa::Edge> &edges,
                                           std::vector<HypernodeID> &to_hypergraph,
                                           std::vector<HypernodeID> &to_graph) {

  for (kaffpa::NodeID u = 0; u + 1 < nodes.size(); ++u) {
    const auto &node = nodes[u];
    const auto &next = nodes[u + 1];
    const auto weight = node.weight;
    const auto num_predecessors = next.firstInEdge - node.firstInEdge;
    if (weight > 0) {
      to_graph.push_back(u);
    } else if (num_predecessors == 0) {
      to_graph.push_back(u);
    }
  }

  to_hypergraph.resize(nodes.size() - 1, Hypergraph::kInvalidHypernodeID);
  for (std::size_t i = 0; i < to_graph.size(); ++i) {
    to_hypergraph[to_graph[i]] = i;
  }
}

static inline void updatePartitionTableAndResult(const std::string &filename,
                                                 const kaffpa::KaffpaHeader &header,
                                                 const std::vector<kaffpa::PartitionID> &table,
                                                 const kaffpa::KaffpaResult &result) {
  FILE *shm = std::fopen(filename.c_str(), "r+");

  std::size_t pos = sizeof(kaffpa::KaffpaHeader)
                    + (header.numberOfNodes + 1) * sizeof(kaffpa::Node)
                    + header.numberOfEdges * sizeof(kaffpa::Edge)
                    + header.numberOfEdges * sizeof(kaffpa::InEdge)
                    + sizeof(kaffpa::PartitionConfig);

  std::fseek(shm, pos, SEEK_SET);

  std::fwrite(table.data(), sizeof(kaffpa::PartitionID), header.numberOfNodes, shm);
  std::fwrite(&result, sizeof(kaffpa::KaffpaResult), 1, shm);
  std::fclose(shm);
}

static inline void writePartitionToSharedMemoryGraphFile(const Hypergraph &hg, const std::string &filename) {
  // compute edge cut
  kaffpa::EdgeWeight cut = 0;
  for (const HypernodeID &u : hg.nodes()) {
    for (const HyperedgeID &e : hg.incidentHeadEdges(u)) {
      for (const HypernodeID &v : hg.tails(e)) {
        if (hg.partID(u) != hg.partID(v)) {
          cut += hg.edgeWeight(e);
        }
      }
    }
  }

  kaffpa::KaffpaResult result{
      .k =  hg.k(),
      .edgeCut = cut,
      .objective = cut,
      .errorCount = 0,
      .warningCount = 0,
      .numRounds = 1
  };

  // build partition table for the original graph
  kaffpa::KaffpaHeader header{};
  std::vector<kaffpa::Node> nodes;
  std::vector<kaffpa::Edge> forward_edges;
  std::vector<kaffpa::InEdge> backward_edges;
  kaffpa::PartitionConfig config{};
  std::vector<HypernodeID> to_hypergraph;
  std::vector<HypernodeID> to_graph;

  kahypar::io::readBinaryKaffpaD(filename, header, nodes, forward_edges, backward_edges, config);
  kahypar::io::buildBinaryKaffpaDHNMap(nodes, forward_edges, to_hypergraph, to_graph);

  std::vector<kaffpa::PartitionID> table;
  bool handled_input_out = false;

  for (kaffpa::NodeID u = 0; u < header.numberOfNodes; ++u) {
    ASSERT(u + 1 < nodes.size());



    const auto &node = nodes[u];
    const auto &next = nodes[u + 1];
    const auto weight = node.weight;
    const auto num_predecessors = next.firstInEdge - node.firstInEdge;

    if (node.weight2 != 0) {
      throw std::runtime_error("bad kaffpaD format: weight2 is not null for some node");
    }

    if (weight > 0) {
      ASSERT(to_hypergraph[u] != Hypergraph::kInvalidHypernodeID);
      table.push_back(hg.partID(to_hypergraph[u]));
    } else if (!handled_input_out && num_predecessors == 0) {
      ASSERT(to_hypergraph[u] != Hypergraph::kInvalidHypernodeID);
      table.push_back(hg.partID(to_hypergraph[u]));
      handled_input_out = true;
    } else {
      ASSERT(num_predecessors == 1);
      const auto pred_id = backward_edges[node.firstInEdge].target;
      ASSERT(nodes[pred_id].weight + nodes[pred_id].weight2 > 0);
      ASSERT(to_hypergraph[pred_id] != Hypergraph::kInvalidHypernodeID);
      table.push_back(hg.partID(to_hypergraph[pred_id]));
    }
  }

  // validate edge cut TODO
//  ASSERT([&]() {
//    kaffpa::EdgeWeight _validate_cut = 0;
//    for (kaffpa::NodeID u = 0; u < header.numberOfNodes; ++u) {
//      const auto &node = nodes[u];
//      const auto &next = nodes[u + 1];
//      for (kaffpa::EdgeID e = node.firstOutEdge; e < next.firstOutEdge; ++e) {
//        kaffpa::NodeID v = forward_edges[e].target;
//        if (table[u] != table[v]) {
//          _validate_cut += forward_edges[e].weight;
//        }
//      }
//    }
//    return _validate_cut == cut;
//  }());

  // permute block ids so that they're in topological order
  std::vector<std::vector<bool>> Q(result.k, std::vector<bool>(result.k));
  for (kaffpa::NodeID u = 0; u < header.numberOfNodes; ++u) {
    const auto &node = nodes[u];
    const auto &next = nodes[u + 1];
    const auto part_u = table[u];

    for (kaffpa::EdgeID e = node.firstOutEdge; e < next.firstOutEdge; ++e) {
      const kaffpa::NodeID v = forward_edges[e].target;
      const auto part_v = table[v];
      if (part_u != part_v) {
        Q[part_u][part_v] = true;
      }
    }
  }

  std::vector<kaffpa::NodeID> queue;
  std::vector<std::size_t> permutation(result.k);
  std::vector<std::size_t> in(result.k);
  for (kaffpa::NodeID u = 0; u < result.k; ++u) {
    for (kaffpa::NodeID v = 0; v < result.k; ++v) {
      if (Q[u][v]) {
        ++in[v];
      }
    }
  }
  for (kaffpa::NodeID u = 0; u < result.k; ++u) {
    if (in[u] == 0) {
      queue.push_back(u);
    }
  }
  std::size_t num_popped = 0;
  while (!queue.empty()) {
    const auto u = queue.back();
    queue.pop_back();
    permutation[u] = num_popped;
    ++num_popped;

    for (kaffpa::NodeID v = 0; v < result.k; ++v) {
      if (Q[u][v]) {
        --in[v];
        if (in[v] == 0) {
          queue.push_back(v);
        }
      }
    }
  }
  if (num_popped != result.k) {
    throw std::runtime_error("cannot permute block ids into topological order because the partition is cyclic!");
  }

  for (kaffpa::NodeID u = 0; u < header.numberOfNodes; ++u) {
    table[u] = permutation[table[u]];
  }

  kahypar::io::updatePartitionTableAndResult(filename, header, table, result);
}

static inline Hypergraph createHypergraphFromSharedMemoryGraphFile(const std::string &filename,
                                                                   const PartitionID num_parts) {
  kaffpa::KaffpaHeader header{};
  std::vector<kaffpa::Node> nodes;
  std::vector<kaffpa::Edge> forward_edges;
  std::vector<kaffpa::InEdge> backward_edges;
  kaffpa::PartitionConfig config{};
  kahypar::io::readBinaryKaffpaD(filename, header, nodes, forward_edges, backward_edges, config);

  std::vector<HypernodeID> to_hypergraph;
  std::vector<HypernodeID> to_graph;
  buildBinaryKaffpaDHNMap(nodes, forward_edges, to_hypergraph, to_graph);

  HypernodeID num_hypernodes = 0;
  HyperedgeID num_hyperedges = 0;
  for (kaffpa::NodeID u = 0; u < header.numberOfNodes; ++u) {
    const auto &node = nodes[u];
    const auto &next = nodes[u + 1];
    const auto weight = node.weight;
    const auto num_successors = next.firstOutEdge - node.firstOutEdge;
    const auto num_predecessors = next.firstInEdge - node.firstInEdge;

    if (weight > 0) {
      ++num_hypernodes;
      num_hyperedges += num_successors;
    } else if (num_predecessors == 0) { // input node
      ++num_hypernodes;
      ++num_hyperedges;
    }
  }

  HyperedgeIndexVector index_vector;
  HyperedgeVector edge_vector;
  HypernodeWeightVector hypernode_weights;
  HyperedgeWeightVector hyperedge_weights;
  NumHeadsVector num_heads_vector;

  index_vector.push_back(edge_vector.size());

  for (kaffpa::NodeID u = 0; u < header.numberOfNodes; ++u) {
    const auto &node = nodes[u];
    const auto &next = nodes[u + 1];
    const auto weight = node.weight;
    const auto num_successors = next.firstOutEdge - node.firstOutEdge;
    const auto num_predecessors = next.firstInEdge - node.firstInEdge;
    if (weight == 0 && num_predecessors > 0) {
      continue;
    }
    if (weight == 0 && num_predecessors == 0) {
      edge_vector.push_back(to_hypergraph[u]);
      const auto &edge = forward_edges[node.firstOutEdge];

      for (kaffpa::NodeID e = node.firstOutEdge; e < next.firstOutEdge; ++e) {
        ASSERT(forward_edges[e].weight == edge.weight);
        const auto v = forward_edges[e].target;
        ASSERT(to_hypergraph[v] != Hypergraph::kInvalidHypernodeID);
        edge_vector.push_back(to_hypergraph[v]);
      }

      hyperedge_weights.push_back(edge.weight);
      num_heads_vector.push_back(1);
      index_vector.push_back(edge_vector.size());
    } else if (weight > 0 && num_successors > 0) {
      for (kaffpa::EdgeID aux_e = node.firstOutEdge; aux_e < next.firstOutEdge; ++aux_e) {
        const auto &edge = forward_edges[aux_e];
        const auto &aux_node = nodes[edge.target];
        const auto &next_aux_node = nodes[edge.target + 1];
        ASSERT(aux_node.weight + aux_node.weight2 == 0);

        hyperedge_weights.push_back(edge.weight);

        ASSERT(to_hypergraph[u] != Hypergraph::kInvalidHypernodeID);
        edge_vector.push_back(to_hypergraph[u]);

        for (kaffpa::NodeID e = aux_node.firstOutEdge; e < next_aux_node.firstOutEdge; ++e) {
//          ASSERT(forward_edges[e].weight == edge.weight); // TODO
          const auto v = forward_edges[e].target;
          ASSERT(to_hypergraph[v] != Hypergraph::kInvalidHypernodeID);
          edge_vector.push_back(to_hypergraph[v]);
        }

        num_heads_vector.push_back(1);
        index_vector.push_back(edge_vector.size());
      }
    }

    hypernode_weights.push_back(weight);
  }

  const bool is_directed = true;

  return Hypergraph(num_hypernodes, num_hyperedges, index_vector, edge_vector,
                    is_directed, num_heads_vector, num_parts,
                    &hyperedge_weights, &hypernode_weights);
}
}  // namespace io
}  // namespace kahypar
