#include <fstream>
#include <iostream>
#include <string>

#include "kahypar/definitions.h"
#include "kahypar/io/hypergraph_io.h"
#include "kahypar/utils/string.h"

#include "hypergraph_checker.h"

using namespace std::string_literals;
using namespace kahypar;

std::pair<std::string, std::vector<std::string>> parse_line(const std::string& line) {
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

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cout << "No .bench.txt file specified" << std::endl;
    std::exit(0);
  }

  std::string iscas_filename = argv[1];
  std::string hgr_filename = iscas_filename + ".hgr";

  std::string line;
  std::ifstream iscas(iscas_filename);
  std::unordered_map<std::string, std::vector<std::string>> graph;
  while (std::getline(iscas, line)) {
    if (line[0] == '#') {
      continue;
    }

    auto head_and_tails = parse_line(line);
    auto& head = head_and_tails.first;
    auto& tails = head_and_tails.second;
    if (!head.empty()) {
      graph[std::move(head)] = std::move(tails);
    }
  }
  iscas.close();

  std::unordered_map<std::string, std::size_t> name_to_id;
  std::size_t next_id = 1; // ids start at 1 in the hypergraph file format
  std::size_t num_hypernodes = graph.size();
  std::size_t num_hyperedges = 0;
  for (const auto& pair : graph) {
    name_to_id[pair.first] = next_id++;
    if (!pair.second.empty()) ++num_hyperedges;
  }

  // write hgr
  std::ofstream out(hgr_filename);
  std::vector<std::size_t> pins;
  out << num_hyperedges << " " << num_hypernodes << "\n";
  for (const auto& pair : graph) {
    if (pair.second.empty()) continue; // no tails
    std::size_t pin = name_to_id[pair.first];
    out << pin << " "; // head
    pins.push_back(pin);

    for (const auto& tail : pair.second) { // tails
      pin = name_to_id[tail];
      // TODO is there a better approach to model multi-edge instead of just ignoring them?
      if (std::find(pins.begin(), pins.end(), pin) == pins.end()) { // prevent multi-pins that result from multi-edges
        out << pin << " ";
        pins.push_back(pin);
      }
    }

    out << "\n";
    pins.clear();
  }
  out.close();

  // check result
  try {
    validateHypergraphFile(hgr_filename);
  } catch (const BadHypernodeException& e) {
    std::cout << "Error: " << e.what() << std::endl;
    const HypernodeID& hn = e.hypernode();
    for (const auto& pair : name_to_id) {
      if (pair.second == hn) {
        std::cout << "Note: " << pair.first << " is mapped to " << hn << std::endl;
        break;
      }
    }
    return -1;
  } catch (const BadHypergraphException& e) {
    std::cout << "Error: " << e.what() << std::endl;
    return -1;
  }
  return 0;
}