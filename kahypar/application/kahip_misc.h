#pragma once

namespace kahip {
// gathered from kaffpawrapper.cpp, must be adjusted

struct KaffpaHeader {
  uint64_t numberOfNodes;
  uint64_t numberOfEdges;
  bool hasSecondNodeWeight;
  bool hasNodeFitness;
};

struct KaffpaResult {
  int64_t edgeCut;
  int errorCount;
  int warningCount;
};

// copied from graph_access.h

struct Node {
  EdgeID firstOutEdge;
  InEdgeID firstInEdge;
  NodeWeight weight;
};

struct Edge {
  NodeID target;
  EdgeWeight weight;
};

struct InEdge {
  NodeID target;
  EdgeID id;
};
}