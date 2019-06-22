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
  std::string current_head;
  for (std::size_t pin = 0; pin < num_pins; ++pin) {
    std::string module, type, direction;
    in >> module >> type >> direction;
    if (type == "s") {
      current_head = module;
    } else {
      ASSERT(type == "l", V(pin) << V(module) << V(type) << V(direction));
      ASSERT(!current_head.empty(), V(pin) << V(module) << V(type) << V(direction));
      auto& tails = map[current_head];
      if (std::find(tails.begin(), tails.end(), module) == tails.end()) {
        tails.push_back(module);
      }
    }
  }
  in.close();

  for (const auto& entry : map) {
    LOG << entry.first << " -- " << entry.second;
  }

  std::unordered_map<std::string, std::size_t> name_to_id;
  std::size_t next_id = 1;
  for (const auto& entry : map) {
    const std::string& name = entry.first;
    if (name_to_id[name] == 0) {
      LOG << name << "->" << next_id;
      name_to_id[name] = next_id++;
    }
    for (const std::string& tail : entry.second) {
      if (name_to_id[tail] == 0) {
        LOG << tail << "->" << next_id;
        name_to_id[tail] = next_id++;
      }
    }
  }

  std::ofstream out(hgr_filename);
  const std::size_t num_hyperedges = map.size();
  const std::size_t num_hypernodes = next_id;
  out << num_hyperedges << " " << num_hypernodes << " 11 1\n";
  for (const auto& entry : map) {
    const std::string& head = entry.first;
    const std::vector<std::string>& tails = entry.second;
    out << "1 1 " << name_to_id[head] << " ";
    for (const std::string& tail : tails) {
      out << name_to_id[tail] << " ";
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

