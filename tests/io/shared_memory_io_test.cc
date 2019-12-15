#include <sys/mman.h>
#include "gmock/gmock.h"

#include "kahypar/io/hypergraph_io.h"

void run_test(const std::string& name) {
  std::string path = "/Users/danielseemaier/HyperDAG/tmp/local_laplacian_";
  path += name;
  path += ".graph";

  kahypar::Hypergraph hg(
      kahypar::io::createHypergraphFromSharedMemoryGraphFile(path, 8)
  );
  for (const kahypar::HypernodeID &hn : hg.nodes()) {
    hg.setNodePart(hn, hn % 4);
  }
  kahypar::io::writePartitionToSharedMemoryGraphFile(hg, path);
}

TEST(SharedMemoryIO, LoadGraph) {
  for (int i = 4; i <= 10; ++i) {
    for (int j = 4; j <= 10; ++j) {
      run_test(std::to_string(i) + "_" + std::to_string(j));
    }
  }
}

