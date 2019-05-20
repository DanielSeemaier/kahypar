#pragma once

#include <algorithm>
#include <numeric>
#include <vector>

#include "kahypar/macros.h"

namespace kahypar {
namespace ds {
/**!
 * A dynamic ordering stores an order on some number of elements that can be queried and updated. More precisely, the
 * data structure can answer queries "x < y" and "x > y" and can move objects before or after another object.
 */
class DynamicOrdering {
 public:
  /**!
   * Answers "x < y?"
   *
   * \param x First object
   * \param y Second object
   * \return Whether x < y.
   */
  virtual bool before(std::size_t x, std::size_t y) = 0;

  /**!
   * Answers "x > y?"
   *
   * \param x First object
   * \param y Second object
   * \return Whether x > y.
   */
  bool after(std::size_t x, std::size_t y) {
    return !before(x, y);
  }

  /**!
   * Moves an object just before another object.
   *
   * \param x Object to be moved
   * \param y Reference object
   */
  virtual void moveBefore(std::size_t x, std::size_t y) = 0;

  /**!
   * Moves an object just after another object.
   *
   * \param x Object to be moved
   * \param y Reference object
   */
  virtual void moveAfter(std::size_t x, std::size_t y) = 0;
};

/**!
 * Naive implementation that uses consecutive positions stored in a vector. This representation allows constant time
 * queries but needs linear time for updates.
 */
class NaiveDynamicOrdering : public DynamicOrdering {
 public:
  explicit NaiveDynamicOrdering(std::size_t n) :
    _topological_ordering(n) {
    std::iota(_topological_ordering.begin(), _topological_ordering.end(), 0);
  }

  NaiveDynamicOrdering(const NaiveDynamicOrdering& other) = delete;
  NaiveDynamicOrdering& operator=(const NaiveDynamicOrdering& other) = delete;

  NaiveDynamicOrdering(NaiveDynamicOrdering&& other) = default;
  NaiveDynamicOrdering& operator=(NaiveDynamicOrdering&& other) = default;

  ~NaiveDynamicOrdering() = default;

  bool before(const std::size_t x, const std::size_t y) override {
    ASSERT(x < _topological_ordering.size());
    ASSERT(y < _topological_ordering.size());
    return _topological_ordering[x] < _topological_ordering[y];
  }

  void moveBefore(std::size_t x, std::size_t y) override {
    const std::size_t pos_x = _topological_ordering[x];
    const std::size_t pos_y = _topological_ordering[y];

    if (before(x, y)) {
      for (std::size_t& position : _topological_ordering) {
        if (pos_x < position && position < pos_y) {
          --position;
        }
      }
    } else {
      for (std::size_t& position : _topological_ordering) {
        if (pos_y <= position && position < pos_x) {
          ++position;
        }
      }
    }
    ASSERT(_topological_ordering[y] > 0);
    _topological_ordering[x] = _topological_ordering[y] - 1;
  }

  void moveAfter(std::size_t x, std::size_t y) override {
    moveBefore(x, y);
    std::swap(_topological_ordering[x], _topological_ordering[y]);
  }

 private:
  std::vector<std::size_t> _topological_ordering;
};
} // namespace ds

using ds::DynamicOrdering;
using ds::NaiveDynamicOrdering;
} // namespace kahypar