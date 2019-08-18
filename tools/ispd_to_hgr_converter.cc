#include "hypergraph_checker.h"

using namespace kahypar;

int main(int argc, char* argv[]) {
  if (argc != 3) {
    std::cout << "Usage: " << argv[0] << " <*.netD> <*.hgr>" << std::endl;
    std::exit(0);
  }

  std::string netd_filename = argv[1];
  std::string hgr_filename = argv[2];

  std::ifstream in(netd_filename);
  std::size_t ignored, num_pins, num_nets, num_modules, pad_offset;
  in >> ignored >> num_pins >> num_nets >> num_modules >> pad_offset;

  std::unordered_map<std::string, std::vector<std::string>> map;
  std::vector<std::pair<std::vector<std::string>, std::vector<std::string>>> graph;
  std::string current_head;
  for (std::size_t pin = 0; pin < num_pins; ++pin) {
    std::string module, type, direction;
    in >> module >> type >> direction;
    ASSERT(type == "s" || type == "l");
    ASSERT(module[0] == 'a' || module[0] == 'p');
    ASSERT(direction == "B" || direction == "I" || direction == "O");

    if (type == "s") {
      if (!graph.empty()) {
        if (graph.back().first.size() + graph.back().second.size() <= 1) {
          graph.pop_back();
        }
      }
      graph.emplace_back();
    }

    ASSERT(!graph.empty());

    // skip pads and bidirectional signals
    if (module[0] == 'p' || direction[0] == 'B') {
      continue;
    }

    if (direction == "I") {
      graph.back().first.push_back(module);
    } else {
      ASSERT(direction == "O");
      graph.back().second.push_back(module);
    }
  }
  in.close();

  std::unordered_map<std::string, std::size_t> name_to_id;
  std::size_t next_id = 1;
  for (const auto& edge : graph) {
    for (const auto& head : edge.first) {
      if (name_to_id[head] == 0) {
        LOG << head << "->" << next_id;
        name_to_id[head] = next_id++;
      }
    }
    for (const auto& tail : edge.second) {
      if (name_to_id[tail] == 0) {
        LOG << tail << "->" << next_id;
        name_to_id[tail] = next_id++;
      }
    }
  }

  auto name_to_id_f = [](const std::string &name) {
    return std::strtol(name.data() + 1, nullptr, 10) + 1;
  };

  std::ofstream out(hgr_filename);
  const std::size_t num_hyperedges = graph.size();
  const std::size_t num_hypernodes = next_id - 1;
  out << num_hyperedges << " " << num_hypernodes << " 11 1\n";
  for (const auto& edge : graph) {
    out << "1 "; // hyperedge weight
    out << edge.first.size() << " "; // num heads
    for (const auto& head : edge.first) {
      out << name_to_id_f(head) << " ";
    }
    for (const auto& tail : edge.second) {
      out << name_to_id_f(tail) << " ";
    }
    out << "\n";
  }

  for (std::size_t i = 0; i < num_hypernodes; ++i) {
    out << "1\n";
  }
  out.close();

  validateHypergraphFile(hgr_filename);
  return 0;
}

