#pragma once

#include "kahypar/partition/dhgp/topord.h"

#ifndef KAHYPAR_ENABLE_DHGP
#error This file should not be included when building without KAHYPAR_ENABLE_DHGP
#endif // KAHYPAR_ENABLE_DHGP

namespace kahypar::dhgp {
static inline void partition(Hypergraph& hypergraph, const Context& context) {

}
}