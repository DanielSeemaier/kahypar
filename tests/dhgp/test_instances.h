#pragma once

#include "kahypar/datastructure/hypergraph.h"

namespace kahypar::dhgp::test_instances {
Hypergraph c17() {
  return Hypergraph(11, 6, HyperedgeIndexVector { 0, 3, 6, 9, 12, 15, 18 },
                    HyperedgeVector { 0, 2, 7, 1, 8, 2, 2, 10, 4, 3, 5, 1, 5, 6, 10, 9, 1, 0 },
                    true, HyperedgeVector { 1, 1, 1, 1, 1, 1 });
}
}