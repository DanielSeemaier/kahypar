#include "hypergraph_checker.h"

#include <kahypar/utils/stringhelpers.h>

using namespace kahypar;

int main(int argc, char* argv[]) {
  if (argc != 3) {
    std::cout << "Usage: " << argv[0] << " <*.nets> <*.hgr>" << std::endl;
    std::exit(0);
  }

  std::string nets_filename = argv[1];
  std::string hgr_filename = argv[2];

  std::ifstream in(nets_filename);
  std::string line;

  LOG << "Reading *.nets file ...";
  bool read_net = false;
  std::string net_name;
  std::vector<std::vector<std::pair<std::string, bool>>> graph;
  while (std::getline(in, line)) {
    line = string::trim(line);
    if (string::starts_with(line, "NetDegree")) {
      read_net = true;
      std::string tmp;
      std::stringstream(line) >> tmp >> tmp >> tmp >> net_name;
      graph.emplace_back();
    } else if (read_net) {
      std::string pin, type;
      std::stringstream(line) >> pin >> type;
      if (graph.empty()) {
        LOG << "Bad file format";
        std::exit(1);
      }
      bool add_edge = true;
      for (const auto& prev_edge : graph.back()) {
        if (prev_edge.first == pin) {
//          LOG << "Warning: edge with multi-pins:" << net_name;
//          for (const auto& prev_edge_prime : graph.back()) {
//            LOG << "\t" << prev_edge_prime.first << prev_edge_prime.second;
//          }
//          LOG << "\t" << pin << (type == "O");
          add_edge = false;
          break;
        }
      }
      if (add_edge) {
        graph.back().emplace_back(pin, type == "O");
      }
    }
  }
  in.close();

  LOG << "Mapping node IDs to integers ...";
  std::unordered_map<std::string, std::size_t> name_to_id;
  std::size_t next_id = 1;
  for (const auto& edge : graph) {
    for (const auto& pin : edge) {
      const auto& name = pin.first;
      if (name_to_id[name] == 0) {
        name_to_id[name] = next_id++;
      }
    }
  }

  LOG << "Writing *.hgr file ...";
  std::ofstream out(hgr_filename);
  out << graph.size() << " " << next_id - 1 << " 11 1\n";

  std::size_t total_num_pins = 0;
  std::size_t total_num_heads = 0;
  std::size_t total_num_tails = 0;
  for (const auto& edge : graph) {
    std::size_t num_heads = 0;
    for (const auto& pin : edge) {
      ++total_num_pins;
      if (pin.second) {
        ++total_num_heads;
        ++num_heads;
      } else {
        ++total_num_tails;
      }
    }

    out << "1 " << num_heads << " ";
    for (const auto& pin : edge) {
      if (pin.second) {
        out << name_to_id[pin.first] << " ";
      }
    }

    for (const auto& pin : edge) {
      if (!pin.second) {
        out << name_to_id[pin.first] << " ";
      }
    }

    out << "\n";
  }

  for (std::size_t i = 0; i < next_id; ++i) {
    out << "1\n";
  }
  out.close();

  LOG << "# HEs:" << graph.size();
  LOG << "# Pins:" << total_num_pins;
  LOG << "# Heads:" << total_num_heads;
  LOG << "# Tails:" << total_num_tails;

  validateHypergraphFile(hgr_filename);
  return 0;
}

