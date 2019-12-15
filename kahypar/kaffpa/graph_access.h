/******************************************************************************
 * graph_access.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#pragma once

#include <vector>

#include "kahypar/kaffpa/definitions.h"

namespace kaffpa {

using QString = std::string;

struct Node {
  EdgeID firstOutEdge;
  InEdgeID firstInEdge;
  NodeWeight weight;
  NodeWeight weight2;
  int fitness;
};

struct Edge {
  NodeID target;
  EdgeWeight weight;
  bool isForcedCut;
};

struct InEdge {
  NodeID target;
  EdgeID id;
};

class KaffpaWrapperInternal;

class BasicGraph {
  friend class KaffpaWrapperInternal;

public:
  BasicGraph(QString name);

  virtual ~BasicGraph();

  /* graph access methods */
  NodeID getNumNodes() const { return mNodes.size() - 1; }

  EdgeID getNumEdges() const { return mEdgesForward.size(); }

  EdgeID get_first_out_edge(NodeID node) const { return mNodes[node].firstOutEdge; }

  EdgeID get_first_invalid_out_edge(NodeID node) const { return mNodes[node + 1].firstOutEdge; }

  InEdgeID get_first_in_edge(NodeID node) const { return mNodes[node].firstInEdge; }

  InEdgeID get_first_invalid_in_edge(NodeID node) const { return mNodes[node + 1].firstInEdge; }

  NodeWeight getNodeWeight(NodeID node) const { return mNodes[node].weight; }

  void setNodeWeight(NodeID node, NodeWeight weight) { mNodes[node].weight = weight; }

  bool hasNonDefaultNodeWeight() const { return mHasNonDefaultNodeWeight; }

  void setHasNonDefaultNodeWeight(bool b) { mHasNonDefaultNodeWeight = b; }

  NodeWeight getNodeWeight2(NodeID node) const { return mNodes[node].weight2; }

  void setNodeWeight2(NodeID node, NodeWeight weight) { mNodes[node].weight2 = weight; }

  bool hasSecondNodeWeight() const { return mHasSecondNodeWeight; }

  void setHasSecondNodeWeight(bool b) { mHasSecondNodeWeight = b; }

  int getNodeFitness(NodeID node) const { return mNodes[node].fitness; }

  void setNodeFitness(NodeID node, int val) { mNodes[node].fitness = val; }

  bool hasNodeFitness() const { return mHasNodeFitness; }

  void setHasNodeFitness(bool b) { mHasNodeFitness = b; }

  EdgeID getNodeDegree(NodeID node) const;

  EdgeID getNodeInDegree(NodeID node) const;

  EdgeWeight getWeightedNodeDegree(NodeID node) const;

  EdgeWeight getEdgeWeight(EdgeID edge) const { return mEdgesForward[edge].weight; }

  void setEdgeWeight(EdgeID edge, EdgeWeight weight) { mEdgesForward[edge].weight = weight; }

  EdgeWeight getInEdgeWeight(InEdgeID edge) const { return mEdgesForward[mEdgesBackwards[edge].id].weight; }

  bool hasNonDefaultEdgeWeight() const { return mHasNonDefaultEdgeWeight; }

  void setHasNonDefaultEdgeWeight(bool b) { mHasNonDefaultEdgeWeight = b; }

  NodeID getEdgeTarget(EdgeID edge) const { return mEdgesForward[edge].target; }

  NodeID getEdgeSource(InEdgeID edge) const { return mEdgesBackwards[edge].target; }

  NodeID edgeIsForcedCut(EdgeID edge) const { return mEdgesForward[edge].isForcedCut; }

  void setEdgeForcedCut(EdgeID edge) { mEdgesForward[edge].isForcedCut = true; }

  NodeID inEdgeIsForcedCut(InEdgeID edge) const { return mEdgesForward[mEdgesBackwards[edge].id].isForcedCut; }

  EdgeRatingType getEdgeRating(EdgeID edge) const { return mCoarseningEdgeProps[edge]; }

  void setEdgeRating(EdgeID edge, EdgeRatingType rating) { mCoarseningEdgeProps[edge] = rating; }

  /* partition indices */
  PartitionID getNumPartitions() const { return mNumPartitions; }

  void setNumPartitions(PartitionID count) { mNumPartitions = count; }

  PartitionID getPartitionOfNode(NodeID node) const { return mPartitionTable[node]; }

  void setPartitionOfNode(NodeID node, PartitionID id) { mPartitionTable[node] = id; }

  PartitionID getSecondPartitionOfNode(NodeID node) const { return mSecondPartitionTable[node]; }

  void setSecondPartitionOfNode(NodeID node, PartitionID id) { mSecondPartitionTable[node] = id; }

  PartitionID getOriginalPartitionOfNode(NodeID node) const { return mOriginalPartitionTable[node]; }

  void setOriginalPartitionOfNode(NodeID node, PartitionID id) { mOriginalPartitionTable[node] = id; }

  void resizeSecondPartitionTable(int n) { mSecondPartitionTable.resize(n); }

  void resizeOriginalPartitionTable(int n) { mOriginalPartitionTable.resize(n); }

  /* build methods */
  void startConstruction(NodeID nodes, EdgeID edges);

  NodeID newNode();

  EdgeID newEdge(NodeID source, NodeID target);

  bool finishConstruction();

  NodeID getCurrentNode() const { return mNodeIdx; }

  EdgeID getCurrentEdge() const { return mEdgeIdx; }

  void dump_node(NodeID y);

  QString print(int mode = 0) const { return printGraph(0, mode == 1); }

  QString printGraph(int useList = 0, bool printOriginalPartitioning = false) const;

private:
  QString mName;

  NodeID mNumNodes;
  NodeID mNumEdges;
  std::vector<Node> mNodes;
  std::vector<Edge> mEdgesForward;
  std::vector<InEdge> mEdgesBackwards;
  std::vector<EdgeRatingType> mCoarseningEdgeProps;
  bool mHasNonDefaultNodeWeight;
  bool mHasSecondNodeWeight;
  bool mHasNodeFitness;
  bool mHasNonDefaultEdgeWeight;

  /* construction properties */
  bool mBuildingGraph;
  int mLastSource;
  NodeID mNodeIdx; /* current node that is constructed */
  EdgeID mEdgeIdx; /* current edge that is constructed */

  PartitionID mNumPartitions;
  std::vector<PartitionID> mPartitionTable;
  std::vector<PartitionID> mSecondPartitionTable;
  std::vector<PartitionID> mOriginalPartitionTable;
};

} // namespace kaffpa