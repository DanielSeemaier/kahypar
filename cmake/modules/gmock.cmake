# taken from http://johnlamp.net/cmake-tutorial-5-functionally-improved-testing.html
function(add_gmock_test target)
    add_executable(${target} ${ARGN})
    target_link_libraries(${target} gtest gmock_main ${CMAKE_THREAD_LIBS_INIT})

    set_property(TARGET ${target} PROPERTY CXX_STANDARD 14)
    set_property(TARGET ${target} PROPERTY CXX_STANDARD_REQUIRED ON)

    add_test(${target} ${target})
endfunction()
