#pragma once

#include "kahypar/datastructure/hypergraph.h"

namespace kahypar::dhgp::test_instances {
Hypergraph c17() {
  return Hypergraph(11, 6, true, HyperedgeIndexVector { 0, 3, 6, 9, 12, 15, 18 },
                    HyperedgeVector { 0, 2, 7, 1, 8, 2, 2, 10, 4, 3, 5, 1, 5, 6, 10, 9, 1, 0 },
                    HyperedgeVector { 1, 1, 1, 1, 1, 1 });
}

Hypergraph cycle() {
  return Hypergraph(8, 5, true, HyperedgeIndexVector { 0, 2, 5, 7, 10, 13 },
                    HyperedgeVector { 0, 7, 2, 1, 0, 3, 2, 5, 4, 3, 7, 6, 5 },
                    HyperedgeVector { 1, 1, 1, 1, 1 });
}
}