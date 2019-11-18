#pragma once

#include <stack>

#include "kahypar/definitions.h"

namespace kahypar {

class HTimer {
 public:
  void start() {
    _start_points.push(std::chrono::high_resolution_clock::now());
  }

  double stop() {
    ASSERT(!_start_points.empty());
    const HighResClockTimepoint start = _start_points.top();
    const HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    _start_points.pop();
    return std::chrono::duration<double>(end - start).count();;
  }

  bool running() const {
    return !_start_points.empty();
  }

 private:
  std::stack<HighResClockTimepoint> _start_points;
};

} // namespace kahypar