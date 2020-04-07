#pragma once

#include <algorithm>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "kahypar/definitions.h"

using ::testing::AllOf;
using ::testing::AllOfArray;
using ::testing::Eq;
using ::testing::Matcher;
using ::testing::MatcherInterface;
using ::testing::MatchResultListener;
using ::testing::Property;
using ::testing::UnorderedElementsAreArray;
using ::testing::ResultOf;

namespace kahypar::dhgp::testing {
namespace utility {
template<typename IterPair>
auto iter2vec(const IterPair&& iter_pair) {
  std::vector<typename decltype(iter_pair.first)::value_type> values;
  values.insert(values.begin(), iter_pair.first, iter_pair.second);
  return values;
}
}

auto IsSameDirectedHypergraph(const Hypergraph& as) {
  auto properties_matcher = AllOf(Property(&Hypergraph::initialNumEdges, Eq(as.initialNumEdges())),
                                  Property(&Hypergraph::currentNumEdges, Eq(as.currentNumEdges())),
                                  Property(&Hypergraph::initialNumNodes, Eq(as.initialNumNodes())),
                                  Property(&Hypergraph::currentNumNodes, Eq(as.currentNumNodes())),
                                  Property(&Hypergraph::initialNumPins, Eq(as.initialNumPins())),
                                  Property(&Hypergraph::currentNumPins, Eq(as.currentNumPins())));

  using utility::iter2vec;

  auto create_incident_head_edges_matcher = [&as](const HypernodeID hn) {
    return ResultOf([hn](const Hypergraph& hg) {
      return iter2vec(hg.incidentHeadEdges(hn));
    }, UnorderedElementsAreArray(iter2vec(as.incidentHeadEdges(hn))));
  };
  auto create_incident_tail_edges_matcher = [&as](const HypernodeID hn) {
    return ResultOf([hn](const Hypergraph& hg) {
      return iter2vec(hg.incidentTailEdges(hn));
    }, UnorderedElementsAreArray(iter2vec(as.incidentTailEdges(hn))));
  };

  std::vector<decltype(create_incident_head_edges_matcher(0))> incident_head_edges_matchers;
  std::vector<decltype(create_incident_tail_edges_matcher(0))> incident_tail_edges_matchers;
  for (const HypernodeID hn : as.nodes()) {
    incident_head_edges_matchers.push_back(create_incident_head_edges_matcher(hn));
    incident_tail_edges_matchers.push_back(create_incident_tail_edges_matcher(hn));
  }

  auto create_head_pins_matcher = [&as](const HyperedgeID he) {
    return ResultOf([he](const Hypergraph& hg) {
      return iter2vec(hg.headPins(he));
    }, UnorderedElementsAreArray(iter2vec(as.headPins(he))));
  };
  auto create_tail_pins_matcher = [&as](const HyperedgeID he) {
    return ResultOf([he](const Hypergraph& hg) {
      return iter2vec(hg.tailPins(he));
    }, UnorderedElementsAreArray(iter2vec(as.tailPins(he))));
  };

  std::vector<decltype(create_head_pins_matcher(0))> head_pins_matcher;
  std::vector<decltype(create_tail_pins_matcher(0))> tail_pins_matcher;

  for (const HyperedgeID he : as.edges()) {
    head_pins_matcher.push_back(create_head_pins_matcher(he));
    tail_pins_matcher.push_back(create_tail_pins_matcher(he));
  }

  return AllOf(properties_matcher,
               AllOfArray(incident_head_edges_matchers),
               AllOfArray(incident_tail_edges_matchers),
               AllOfArray(head_pins_matcher),
               AllOfArray(tail_pins_matcher));
}
}