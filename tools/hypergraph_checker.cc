#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <algorithm>

#include "hypergraph_checker.h"

using namespace kahypar;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cout << "Wrong number of arguments!" << std::endl;
    std::cout << "Usage: HypergraphChecker <hypergraph.hgr>" << std::endl;
    return -1;
  }
  std::string filename = argv[1];
  validateHypergraphFileAndPanic(filename);
  return 0;
}
