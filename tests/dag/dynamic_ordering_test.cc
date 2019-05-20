#include "gmock/gmock.h"

#include "kahypar/dag/dynamic_ordering.h"

namespace kahypar {
template<typename Ordering, std::size_t N = 8>
static void _simpleMoveBeforeTest() {
  Ordering ordering(N);

  // 0 < 1 < ... < N
  for (std::size_t i = 0; i < N - 1; ++i) {
    ordering.moveBefore(i, i + 1);
    ASSERT_TRUE(ordering.before(i, i + 1));
  }
  for (std::size_t i = 0; i < N - 1; ++i) {
    ASSERT_TRUE(ordering.before(i, i + 1));
  }
}

template<typename Ordering, std::size_t N = 8>
static void _simpleMoveAfterTest() {
  Ordering ordering(N);

  for (std::size_t i = 0; i < N - 1; ++i) {
    ordering.moveAfter(i, i + 1);
    ASSERT_TRUE(ordering.after(i, i + 1));
  }
}

// NaiveTopologicalOrdering
TEST(NAIVE_DYNAMIC_ORDERING, SimpleMoveBeforeWorks) { _simpleMoveBeforeTest<NaiveDynamicOrdering>(); }
TEST(NAIVE_DYNAMIC_ORDERING, SimpleMoveAfterWorks) { _simpleMoveAfterTest<NaiveDynamicOrdering>(); }
} // namespace kahypar