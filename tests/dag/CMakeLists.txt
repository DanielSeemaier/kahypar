# This test needs test instance files, so we copy them to the corresponding build dir
file(COPY test_instances DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_gmock_test(topological_ordering_test topological_ordering_test.cc)
add_gmock_test(cycle_detector_test cycle_detector_test.cc)
add_gmock_test(quotient_graph_test quotient_graph_test.cc)
add_gmock_test(dynamic_ordering_test dynamic_ordering_test.cc)
add_gmock_test(directed_contraction_test directed_contraction_test.cc)
add_gmock_test(directed_hypergraph_test directed_hypergraph_test.cc)
add_gmock_test(km1_gain_manager_test km1_gain_manager_test.cc)
add_gmock_test(hard_rebalance_test hard_rebalance_test.cc)
add_gmock_test(local_search_test local_search_test.cc)
