#include "gmock/gmock.h"

#include <unordered_set>

#include "kahypar/dag/cycle_detector.h"
#include "kahypar/utils/randomize.h"

using ::testing::Eq;
using ::testing::ContainerEq;

namespace kahypar {
template<typename Detector>
static void _simpleCycleTest() {
  Detector detector(5);
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

template<typename Detector>
static void _anotherSimpleCycleTest() {
  Detector detector(5);

  //   /> 3 <\
  // 0 -> 1 <- 4
  //   \> 2 </
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

template<typename Detector>
static void _rejectSelfLoopTest() {
  Detector detector(2);
  ASSERT_FALSE(detector.connect(0, 0));
  ASSERT_FALSE(detector.connect(1, 1));
}

template<typename Detector, std::size_t N = 512>
static void _bigSimpleCycleTest() {
  Detector detector(N);
  for (std::size_t i = 0; i < N - 1; ++i) {
    ASSERT_TRUE(detector.connect(i, i + 1));
  }
  ASSERT_FALSE(detector.connect(N - 1, 0));
}

template<typename Detector, std::size_t N = 128>
static void _bigTreeTest() {
  Detector detector(N);

  for (std::size_t u = 0; u < N; ++u) {
    for (std::size_t v = u + 1; v < N; ++v) {
      ASSERT_TRUE(detector.connect(u, v));
    }
  }

  for (std::size_t u = N - 2; u > 1; --u) {
    ASSERT_FALSE(detector.connect(N - 1, u));
  }
}

template<typename Detector, std::size_t N = 512, std::size_t M = 1024, int seed = 0>
static void _randomGraphTest() {
  Detector detector(N);
  DFSCycleDetector reference_detector(N);

  auto& random = Randomize::instance();
  random.setSeed(seed);
  std::unordered_set<std::size_t> edges;

  for (std::size_t i = 0; i < M; ++i) {
    std::size_t x = random.getRandomInt(0, N - 1);
    std::size_t y = random.getRandomInt(0, N - 1);

    // don't try to insert the same edge twice
    ASSERT(x < (1u << 15u));
    ASSERT(y < (1u << 15u));
    std::size_t key = (x << 16u) | y;
    if (edges.find(key) != edges.end()) {
      continue;
    }

    bool success_detector = detector.connect(x, y);
    if (success_detector) {
      edges.insert(key);
    }
    bool success_reference_detector = reference_detector.connect(x, y);
    ASSERT_EQ(success_detector, success_reference_detector);
  }
}

// DFSCycle Detector
TEST(DFS_CYCLE_DETECTOR, SimpleCycle) { _simpleCycleTest<DFSCycleDetector>(); }
TEST(DFS_CYCLE_DETECTOR, AnotherSimpleCycle) { _anotherSimpleCycleTest<DFSCycleDetector>(); }
TEST(DFS_CYCLE_DETECTOR, RejectsSelfLoops) { _rejectSelfLoopTest<DFSCycleDetector>(); }
TEST(DFS_CYCLE_DETECTOR, BigSimpleCycle) { _bigSimpleCycleTest<DFSCycleDetector>(); }
TEST(DFS_CYCLE_DETECTOR, BigTree) { _bigTreeTest<DFSCycleDetector>(); }

// KahnCycleDetector
TEST(KAHN_CYCLE_DETECTOR, SimpleCycle) { _simpleCycleTest<KahnCycleDetector>(); }
TEST(KAHN_CYCLE_DETECTOR, AnotherSimpleCycle) { _anotherSimpleCycleTest<KahnCycleDetector>(); }
TEST(KAHN_CYCLE_DETECTOR, RejectsSelfLoops) { _rejectSelfLoopTest<KahnCycleDetector>(); }
TEST(KAHN_CYCLE_DETECTOR, BigSimpleCycle) { _bigSimpleCycleTest<KahnCycleDetector>(); }
TEST(KAHN_CYCLE_DETECTOR, BigTree) { _bigTreeTest<KahnCycleDetector>(); }
TEST(KAHN_CYCLE_DETECTOR, Random) { _randomGraphTest<KahnCycleDetector, 128, 256, 128>(); }
TEST(KAHN_CYCLE_DETECTOR, Random1) { _randomGraphTest<KahnCycleDetector, 256, 512, 256>(); }
TEST(KAHN_CYCLE_DETECTOR, Random2) { _randomGraphTest<KahnCycleDetector, 512, 512, 512>(); }
TEST(KAHN_CYCLE_DETECTOR, Random3) { _randomGraphTest<KahnCycleDetector>(); }

// PseudoTopologicalOrderingCycleDetector
TEST(PSEUDO_TOPO_CYCLE_DETECTOR, SimpleCycle) { _simpleCycleTest<PseudoTopologicalOrderingCycleDetector>(); }
TEST(PSEUDO_TOPO_CYCLE_DETECTOR, AnotherSimpleCycle) { _anotherSimpleCycleTest<PseudoTopologicalOrderingCycleDetector>(); }
TEST(PSEUDO_TOPO_CYCLE_DETECTOR, RejectsSelfLoops) { _rejectSelfLoopTest<PseudoTopologicalOrderingCycleDetector>(); }
TEST(PSEUDO_TOPO_CYCLE_DETECTOR, BigSimpleCycle) { _bigSimpleCycleTest<PseudoTopologicalOrderingCycleDetector>(); }
TEST(PSEUDO_TOPO_CYCLE_DETECTOR, BigTree) { _bigTreeTest<PseudoTopologicalOrderingCycleDetector>(); }
TEST(PSEUDO_TOPO_CYCLE_DETECTOR, Random) { _randomGraphTest<PseudoTopologicalOrderingCycleDetector, 128, 256, 128>(); }
TEST(PSEUDO_TOPO_CYCLE_DETECTOR, Random1) { _randomGraphTest<PseudoTopologicalOrderingCycleDetector, 256, 512, 256>(); }
TEST(PSEUDO_TOPO_CYCLE_DETECTOR, Random2) { _randomGraphTest<PseudoTopologicalOrderingCycleDetector, 512, 512, 512>(); }
TEST(PSEUDO_TOPO_CYCLE_DETECTOR, Random3) { _randomGraphTest<PseudoTopologicalOrderingCycleDetector>(); }
} // namespace kahypar