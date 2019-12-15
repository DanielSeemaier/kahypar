#include <sys/mman.h>
#include <fcntl.h>
#include "gmock/gmock.h"

#include "kahypar/io/hypergraph_io.h"

void readFromShm(const std::string& filename) {
  std::FILE *shm = std::fopen(filename.c_str(), "r");
  if (shm == nullptr) {
    std::perror("Could not open file");
    std::exit(0);
  }

  kaffpa::KaffpaHeader header{};
  std::fread(&header, sizeof (kaffpa::KaffpaHeader), 1, shm);
  LOG << "Nodes:" << header.numberOfNodes << "Edges:" << header.numberOfEdges;

  const std::size_t size_nodes = (header.numberOfNodes + 1) * sizeof (kaffpa::Node);
  const std::size_t size_edges = header.numberOfEdges * sizeof (kaffpa::Edge);
  const std::size_t size_in_edges = header.numberOfEdges * sizeof (kaffpa::InEdge);

  std::vector<kaffpa::Node> nodes(header.numberOfNodes + 1);
  std::vector<kaffpa::Edge> forward_edges(header.numberOfEdges);
  std::vector<kaffpa::InEdge> backward_edges(header.numberOfEdges);
  kaffpa::PartitionConfig config{};
  std::vector<kaffpa::PartitionID> partition_table(header.numberOfNodes);
  kaffpa::KaffpaResult result{};

  const auto num_nodes_read = std::fread(nodes.data(), sizeof (kaffpa::Node), header.numberOfNodes + 1, shm);
  const auto num_forward_edges_read = std::fread(forward_edges.data(), sizeof (kaffpa::Edge), header.numberOfEdges, shm);
  const auto num_backward_edges_read = std::fread(backward_edges.data(), sizeof (kaffpa::InEdge), header.numberOfEdges, shm);
  const auto num_configs_read = std::fread(&config, sizeof (kaffpa::PartitionConfig), 1, shm);
  const auto num_partition_entries_read = std::fread(partition_table.data(), sizeof (kaffpa::PartitionID), header.numberOfNodes, shm);
  const auto num_results_read = std::fread(&result, sizeof (kaffpa::KaffpaResult), 1, shm);

  /**
    memcpy(data, &mHeader, mSizeHeader);
    data += mSizeHeader;

    memcpy(data, mG->mNodes.data(), mSizeNodes);
    data += mSizeNodes;

    memcpy(data, mG->mEdgesForward.data(), mSizeEdges);
    data += mSizeEdges;

    memcpy(data, mG->mEdgesBackwards.data(), mSizeInEdges);
    data += mSizeInEdges;

    memcpy(data, &mConfig, mSizeConfig);
    data += mSizeConfig;

    memcpy(data, mG->mPartitionTable.data(), mSizeTable);
    data += mSizeTable;
   */

  if (num_nodes_read != header.numberOfNodes + 1) {
    LOG << "1";
  }
  if (num_forward_edges_read != header.numberOfEdges) {
    LOG << "2";
  }
  if (num_backward_edges_read != header.numberOfEdges) {
    LOG << "3";
  }
  if (num_configs_read != 1) {
    LOG << "4";
  }
  if (num_partition_entries_read != header.numberOfNodes) {
    LOG << "5";
  }
  if (num_results_read != 1) {
    LOG << "6";
  }
  /*
    mSizeHeader  = sizeof(KaffpaHeader);
    mSizeNodes   = sizeof(Node)        * (mG->getNumNodes() + 1);
    mSizeEdges   = sizeof(Edge)        *  mG->getNumEdges();
    mSizeInEdges = sizeof(InEdge)      *  mG->getNumEdges();
    mSizeConfig  = sizeof(PartitionConfig);
    mSizeTable   = sizeof(PartitionID) *  mG->getNumNodes();
    mSizeResult  = sizeof(KaffpaResult);
   */

  std::fclose(shm);
}

TEST(SharedMemoryIO, LoadGraph) {
  readFromShm("test_instances/test.shm");
}

