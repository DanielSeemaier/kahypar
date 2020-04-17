#include "gmock/gmock.h"

#include "kahypar/partition/dhgp/cycle_detector.h"

using ::testing::TestWithParam;
using ::testing::Values;
using ::testing::Types;

namespace kahypar::dhgp {
template<typename Detector>
class CycleDetectorTest : public ::testing::Test {
};

using cycle_detector_types = Types<
    KahnCycleDetector,
    DFSCycleDetector,
    PseudoTopologicalOrderingCycleDetector>;

TYPED_TEST_CASE(CycleDetectorTest, cycle_detector_types);

TYPED_TEST(CycleDetectorTest, DetectsSimpleCycleAndMaintainsState) {
  TypeParam detector(5);
  // build a path: 0 -> 1 -> 2 -> 3 -> 4
  ASSERT_TRUE(detector.connect(0, 1));
  ASSERT_TRUE(detector.connect(1, 2));
  ASSERT_TRUE(detector.connect(2, 3));
  ASSERT_TRUE(detector.connect(3, 4));

  // these edges should all fail
  ASSERT_FALSE(detector.connect(4, 0));
  ASSERT_FALSE(detector.connect(4, 1));
  ASSERT_FALSE(detector.connect(4, 2));
  ASSERT_FALSE(detector.connect(4, 3));

  // other edges should be ok though
  ASSERT_TRUE(detector.connect(0, 4));
  ASSERT_TRUE(detector.connect(1, 4));
  ASSERT_TRUE(detector.connect(2, 4));
  ASSERT_TRUE(detector.connect(0, 3));
  ASSERT_TRUE(detector.connect(1, 3));
  ASSERT_TRUE(detector.connect(0, 2));
}

TYPED_TEST(CycleDetectorTest, DetectsAnotherCycle) {
  TypeParam detector(5);
  /*   /> 3 <\
     0 -> 1 <- 4
       \> 2 </    */
  ASSERT_TRUE(detector.connect(0, 1));
  ASSERT_TRUE(detector.connect(0, 2));
  ASSERT_TRUE(detector.connect(0, 3));
  ASSERT_TRUE(detector.connect(4, 3));
  ASSERT_TRUE(detector.connect(4, 2));
  ASSERT_TRUE(detector.connect(4, 1));

  // cycles of length 2
  ASSERT_FALSE(detector.connect(3, 4));
  ASSERT_FALSE(detector.connect(2, 4));
  ASSERT_FALSE(detector.connect(1, 4));
  ASSERT_FALSE(detector.connect(1, 0));
  ASSERT_FALSE(detector.connect(2, 0));
  ASSERT_FALSE(detector.connect(3, 0));

  // 4 -> 0, 3 -> 2 is still acyclic though
  ASSERT_TRUE(detector.connect(4, 0));
  ASSERT_TRUE(detector.connect(3, 2));
}

TYPED_TEST(CycleDetectorTest, RejectsSelfLoop) {
  TypeParam detector(5);
  ASSERT_FALSE(detector.connect(0, 0));
  ASSERT_FALSE(detector.connect(1, 1));
}

TYPED_TEST(CycleDetectorTest, DetectsBigCycle) {
  constexpr std::size_t N = 128;
  TypeParam detector(N);

  for (std::size_t i = 0; i < N - 1; ++i) {
    ASSERT_TRUE(detector.connect(i, i + 1));
  }
  ASSERT_FALSE(detector.connect(N - 1, 0));
}

TYPED_TEST(CycleDetectorTest, RejectsCyclesInBigTree) {
  constexpr std::size_t N = 64;
  TypeParam detector(N);

  for (std::size_t u = 0; u < N; ++u) {
    for (std::size_t v = u + 1; v < N; ++v) {
      ASSERT_TRUE(detector.connect(u, v));
    }
  }

  for (std::size_t u = N - 2; u > 1; --u) {
    ASSERT_FALSE(detector.connect(N - 1, u));
  }
}
}// namespace kahypar::dhgp