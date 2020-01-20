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


#include "kahypar/application/command_line_options.h"
#include "kahypar/definitions.h"
#include "kahypar/io/hypergraph_io.h"
#include "kahypar/partitioner_facade.h"
#include "kahypar/partition/context_enum_classes.h"

int main(int argc, char *argv[]) {
  kahypar::Context context;

  kahypar::processCommandLineInput(context, argc, argv);

  context.imbalanced_intermediate_step = false;
  context.reduce_balance_during_uncoarsening = false;

  if (context.bin_kaffpaD) {
    LOG << "Reading input / writing output in kaffpaD shared memory format";
    context.partition.write_partition_file = false;
  }

  kahypar::Hypergraph hypergraph(context.bin_kaffpaD
    ? kahypar::io::createHypergraphFromSharedMemoryGraphFile(context.partition.graph_filename, context.partition.k)
    : kahypar::io::createHypergraphFromFile(context.partition.graph_filename, context.partition.k)
  );

  kahypar::PartitionerFacade().partition(hypergraph, context);

  if (context.partition.mode == kahypar::Mode::acyclic || context.partition.mode == kahypar::Mode::acyclic_kway) {
    if (!kahypar::dag::isAcyclic(hypergraph)) {
      throw std::runtime_error("final partition is cyclic!");
    }
  }

  if (context.bin_kaffpaD) {
    kahypar::io::writePartitionToSharedMemoryGraphFile(hypergraph, context.partition.graph_filename);
  }

  return 0;
}
