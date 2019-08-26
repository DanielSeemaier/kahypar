#include "metis_graph_reader.h"
#include "hypergraph_checker.h"

using namespace kahypar;

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cout << "Usage: " << argv[0] << " <*.graph> <*.hgr> [--skip-verify]" << std::endl;
    std::exit(0);
  }

  std::string graph_filename = argv[1];
  std::string hgr_filename = argv[2];
  const bool skip_verify = (argc == 4) && std::string(argv[3]) == "--skip-verify";

  std::vector<int32_t> node_weights;
  std::vector<int32_t> edge_weights;
  std::vector<std::pair<std::size_t, std::size_t>> edges;
  metis::readMetisGraph(graph_filename, node_weights, edge_weights, edges);

  std::size_t num_hyperedges = 0;
  std::size_t current_head = 0;
  std::size_t num_pins = 0;
  for (const auto& edge : edges) {
    if (edge.first != current_head) {
      if (num_pins > 0) {
        ++num_hyperedges;
      }
      current_head = edge.first;
      num_pins = 0;
    }
    if (edge.first < edge.second) {
      ++num_pins;
    }
  }

  std::ofstream out(hgr_filename);
  out << edge_weights.size() << " " << node_weights.size() << " 11 1\n";

  current_head = 0;
  num_pins = 0;
  bool first_pin = true;
  for (const auto& edge : edges) {
    out << "1 1 " << edge.first + 1 << " " << edge.second + 1 << "\n";
//    if (edge.first != current_head) {
//      if (num_pins > 0) {
//        out << "\n";
//      }
//      current_head = edge.first;
//      num_pins = 0;
//      first_pin = true;
//    }
//    if (edge.first < edge.second) {
//      if (first_pin) {
//        out << "1 1 " << edge.first + 1 << " ";
//        first_pin = false;
//      }
//      out << edge.second + 1 << " ";
//      ++num_pins;
//    }
  }
//  if (num_pins > 0) {
//    out << "\n";
//  }

  for (std::size_t i = 0; i < node_weights.size(); ++i) {
    out << "1\n";
  }
  out.close();

  // check result
  if (!skip_verify) {
    std::cout << "Validating result ... " << std::flush;
    validateHypergraphFileAndPanic(hgr_filename);
    std::cout << "OK\n";
  }
  return 0;
}