function(add_gmock_test_with_gmock target)
  add_gmock_test(${target} ${ARGN})
  target_link_libraries(${target} gmock_main)
endfunction()

add_gmock_test_with_gmock(directed_hypergraph_test directed_hypergraph_test.cc)
add_gmock_test_with_gmock(topord_test topord_test.cc)
add_gmock_test_with_gmock(cycle_detector_test cycle_detector_test.cc)