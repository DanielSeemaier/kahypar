# This test needs test instance files, so we copy them to the corresponding build dir
file(COPY test_instances DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_gmock_test(topological_ordering_test topological_ordering_test.cc)
add_gmock_test(cycle_detector_test cycle_detector_test.cc)