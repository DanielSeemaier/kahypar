#include <sys/mman.h>
#include "gmock/gmock.h"

#include "kahypar/io/hypergraph_io.h"

TEST(SharedMemoryIO, LoadGraph) {
  kahypar::Hypergraph hg(
      kahypar::io::createHypergraphFromSharedMemoryGraphFile("test_instances/test.shm", 4)
  );
  for (const kahypar::HypernodeID &hn : hg.nodes()) {
    hg.setNodePart(hn, hn % 4);
  }

  kahypar::io::writePartitionToSharedMemoryGraphFile(hg, "test_instances/test.shm");
}

