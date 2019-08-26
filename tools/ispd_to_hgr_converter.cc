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
        //LOG << head << "->" << next_id;
        name_to_id[head] = next_id++;
      }
    }
    for (const auto& tail : edge.second) {
      if (name_to_id[tail] == 0) {
        //LOG << tail << "->" << next_id;
        name_to_id[tail] = next_id++;
      }
    }
  }

//  std::unordered_map<std::size_t, std::vector<std::size_t>> other_graph;
//  for (const auto& edge : graph) {
//    for (const auto& head : edge.first) {
//      for (const auto& tail : edge.second) {
//        auto& tails = other_graph[name_to_id[head]];
//        const std::size_t tail_id = name_to_id[tail];
//        if (std::find(tails.begin(), tails.end(), tail_id) == tails.end()) {
//          tails.push_back(tail_id);
//        }
//      }
//    }
//  }

  std::unordered_map<std::size_t, std::pair<std::vector<std::size_t>, std::vector<std::size_t>>> other_graph;
  std::size_t net_id = 0;
  for (const auto& edge : graph) {
    ++net_id;

    for (const auto& head_name : edge.first) {
      const std::size_t head = name_to_id[head_name];
      if (std::find(other_graph[head].second.begin(), other_graph[head].second.end(), net_id) == other_graph[head].second.end()) {
        other_graph[head].second.push_back(net_id);
      }
    }

    for (const auto& tail_name : edge.second) {
      const std::size_t tail = name_to_id[tail_name];
      if (std::find(other_graph[tail].first.begin(), other_graph[tail].first.end(), net_id) == other_graph[tail].first.end()) {
        other_graph[tail].first.push_back(net_id);
      }
    }
  }

  auto name_to_id_f = [](const std::string &name) {
    return std::strtol(name.data() + 1, nullptr, 10) + 1;
  };

  PseudoTopologicalOrderingCycleDetector cd(net_id);
  std::size_t num_rejected_edges = 0;
  std::size_t num_accepted_edges = 0;

  for (const auto& edge : other_graph) {
    const auto& heads = edge.second.first;
    const auto& tails = edge.second.second;
    std::vector<std::size_t> rejected_tails;

    for (const auto& head : heads) {
      for (const auto& tail : tails) {
        if (cd.connect(tail - 1, head - 1)) {
          ++num_accepted_edges;
        } else {
          ++num_rejected_edges;
          if (std::find(rejected_tails.begin(), rejected_tails.end(), tail) == rejected_tails.end()) {
            rejected_tails.push_back(tail);
          }
        }
      }
    }

    auto &nc_tails = other_graph[edge.first].second;
    for (const auto &rejected_tail : rejected_tails) {
      auto pos = std::find(nc_tails.begin(), nc_tails.end(), rejected_tail);
      ASSERT(pos != nc_tails.end(), V(rejected_tails));
      nc_tails.erase(pos);
    }
  }

  LOG << V(num_accepted_edges);
  LOG << V(num_rejected_edges);
  LOG << "fraction_accepted_edges =" << num_accepted_edges / static_cast<double>(num_rejected_edges + num_accepted_edges);

  std::ofstream out(hgr_filename);
//  const std::size_t num_hyperedges = graph.size();
//  const std::size_t num_hypernodes = next_id - 1;
//  out << num_hyperedges << " " << num_hypernodes << " 11 1\n";
//  for (const auto& edge : graph) {
//    out << "1 "; // hyperedge weight
//    out << edge.first.size() << " "; // num heads
//    for (const auto& head : edge.first) {
//      out << name_to_id_f(head) << " ";
//    }
//    for (const auto& tail : edge.second) {
//      out << name_to_id_f(tail) << " ";
//    }
//    out << "\n";
//  }

//  const std::size_t num_hyperedges = other_graph.size();
//  const std::size_t num_hypernodes = next_id - 1;
//  out << num_hyperedges << " " << num_hypernodes << " 11 1\n";
//
//  for (const auto& edge : other_graph) {
//    out << "1 1 " << edge.first << " ";
//    for (const auto& tail : edge.second) {
//      out << tail << " ";
//    }
//    out << "\n";
//  }
//
//  for (std::size_t i = 0; i < num_hypernodes; ++i) {
//    out << "1\n";
//  }
//  out.close();

  std::size_t num_hyperedges = 0;
  for (const auto& edge : other_graph) {
    if (!edge.second.first.empty() && !edge.second.second.empty()) {
      num_hyperedges += 1;//edge.second.first.size();
    }
  }
  const std::size_t num_hypernodes = net_id;
  out << num_hyperedges << " " << num_hypernodes << " 11 1\n";

  for (const auto& edge : other_graph) {
    const auto& heads = edge.second.first;
    const auto& tails = edge.second.second;

    for (const auto& head : heads) {
      out << "1 1 " << head << " ";
      for (const auto& tail : tails) {
        out << tail << " ";
      }
      out << "\n";
      break;
    }
  }

  for (std::size_t i = 0; i < num_hypernodes; ++i) {
    out << "1\n";
  }
  out.close();

  validateHypergraphFile(hgr_filename);
  return 0;
}

