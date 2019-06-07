#include <fstream>
#include <iostream>
#include <string>

#include "kahypar/definitions.h"
#include "kahypar/io/hypergraph_io.h"
#include "kahypar/utils/string.h"

#include "hypergraph_checker.h"

using namespace std::string_literals;
using namespace kahypar;

enum Output {
  HGR,
  GRAPH
};

using Graph = std::unordered_map<std::string, std::vector<std::string>>;
using NameToID = std::unordered_map<std::string, std::size_t>;
using HeadAndTails = std::pair<std::string, std::vector<std::string>>;
using AdjacencyList = std::vector<std::vector<std::size_t>>;

constexpr static bool print_mapping = false;

HeadAndTails parseHeadTails(const std::string& line) {
  // input or output node: INPUT(node) or OUTPUT(node)
  const auto input_length = "INPUT"s.size();
  const auto output_length = "OUTPUT"s.size();
  if (string::starts_with(line, "INPUT")) {
    return {line.substr(input_length + 1, line.size() - input_length - 2), {}};
  }
  if (string::starts_with(line, "OUTPUT")) {
    return {line.substr(output_length + 1, line.size() - output_length - 2), {}};
  }

  // inner node: head = gatter(tail1, tail2, ...)
  auto eq_pos = line.find_first_of('=');
  if (eq_pos == std::string::npos) return {"", {}};

  std::string head = string::trim(line.substr(0, eq_pos));
  auto paren_open_pos = line.find_first_of('(');
  auto paren_close_pos = line.find_first_of(')');
  std::string tails = line.substr(paren_open_pos + 1, paren_close_pos - paren_open_pos - 1);
  return {head, string::split(tails, ",")};
}

void generateHgr(const std::string& out_filename, Graph graph, NameToID name_to_id) {
  const std::size_t num_hypernodes = graph.size();
  std::size_t num_hyperedges = 0;
  for (const auto& pair : graph) {
    if (!pair.second.empty()) {
      ++num_hyperedges;
    }
  }

  std::ofstream out(out_filename);
  if (print_mapping) {
    for (const auto& name_id : name_to_id) {
      out << "% " << name_id.first << "=" << name_id.second << "\n";
    }
  }

  std::vector<std::size_t> pins;
  out << num_hyperedges << " " << num_hypernodes << " 11 1\n";
  for (const auto& pair : graph) {
    if (pair.second.empty()) continue; // no tails
    std::size_t pin = name_to_id[pair.first];
    out << "1 1 " << pin << " "; // head
    pins.push_back(pin);

    for (const auto& tail : pair.second) { // tails
      pin = name_to_id[tail];
      if (std::find(pins.begin(), pins.end(), pin) == pins.end()) { // prevent multi-pins that result from multi-edges
        out << pin << " ";
        pins.push_back(pin);
      }
    }

    out << "\n";
    pins.clear();
  }
  for (std::size_t i = 0; i < num_hypernodes; ++i) {
    out << "1\n";
  }
  out.close();

  // check result
  try {
    validateHypergraphFile(out_filename);
  } catch (const BadHypernodeException& e) {
    std::cout << "Error: " << e.what() << std::endl;
    const HypernodeID& hn = e.hypernode();
    for (const auto& pair : name_to_id) {
      if (pair.second == hn) {
        std::cout << "Note: " << pair.first << " is mapped to " << hn << std::endl;
        break;
      }
    }
    std::exit(-1);
  } catch (const BadHypergraphException& e) {
    std::cout << "Error: " << e.what() << std::endl;
    std::exit(-1);
  }
}

void generateGraph(const std::string& out_filename, Graph graph, NameToID name_to_id) {
  const std::size_t num_nodes = graph.size();
  const std::size_t num_edges = std::accumulate(graph.begin(), graph.end(), 0, [](std::size_t acc, const auto& pair) {
    return acc + pair.second.size();
  });

  std::size_t num_multi_edges = 0;
  AdjacencyList adj_list(num_nodes);
  for (const auto& pair : graph) {
    const auto& head = pair.first;
    const auto& tails = pair.second;
    ASSERT(0 < name_to_id[head] && name_to_id[head] <= num_nodes);
    for (const auto& tail : tails) {
      ASSERT(0 < name_to_id[tail] && name_to_id[tail] <= num_nodes);
      auto& vec = adj_list[name_to_id[tail] - 1];
      auto& to_add = name_to_id[head];
      if (std::find(vec.begin(), vec.end(), to_add) == vec.end()) {
        vec.push_back(to_add);
      } else {
        ++num_multi_edges;
        LOG << "Ignoring multi-edge" << name_to_id[tail] << "-->" << to_add << "===" << tail << "-->" << head;
      }
    }
  }

  std::ofstream out(out_filename);
  if (print_mapping) {
    for (const auto& name_id : name_to_id) {
      out << "% " << name_id.first << "=" << name_id.second << "\n";
    }
  }

  ASSERT(num_edges >= num_multi_edges);
  out << num_nodes << " " << num_edges - num_multi_edges << " 11 1\n";
  for (const auto& heads : adj_list) {
    out << "1 ";
    for (const std::size_t& head : heads) {
      out << head << " 1 ";
    }
    out << "\n";
  }
  out.close();
}

int main(int argc, char* argv[]) {
  if (argc != 3) {
    std::cout << "Usage: " << argv[0] << " <iscas file> <*.hgr|*.graph>" << std::endl;
    std::exit(0);
  }

  std::string iscas_filename = argv[1];
  std::string out_filename = argv[2];
  Output out_format = string::ends_with(out_filename, ".graph")
                      ? Output::GRAPH
                      : Output::HGR;

  std::string line;
  std::ifstream iscas(iscas_filename);
  Graph graph;
  while (std::getline(iscas, line)) {
    if (line[0] == '#') {
      continue;
    }

    auto head_and_tails = parseHeadTails(line);
    auto& head = head_and_tails.first;
    auto& tails = head_and_tails.second;
    if (!head.empty()) {
      graph[std::move(head)] = std::move(tails);
    }
  }
  iscas.close();

  NameToID name_to_id;
  std::size_t next_id = 1; // ids start at 1 in the hypergraph file format
  for (const auto& pair : graph) {
    name_to_id[pair.first] = next_id++;
  }

  switch (out_format) {
    case HGR:
      generateHgr(out_filename, std::move(graph), std::move(name_to_id));
      break;

    case GRAPH:
      generateGraph(out_filename, std::move(graph), std::move(name_to_id));
      break;
  }

  return 0;
}