add_executable(KaHyPar kahypar.cc)
target_link_libraries(KaHyPar ${Boost_LIBRARIES})

add_executable(KaHyParShm kahypar_shm.cc)
target_link_libraries(KaHyParShm ${Boost_LIBRARIES})

add_executable(test2 test.cc)

set_property(TARGET KaHyPar PROPERTY CXX_STANDARD 14)
set_property(TARGET KaHyParShm PROPERTY CXX_STANDARD 14)
set_property(TARGET KaHyPar PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET KaHyParShm PROPERTY CXX_STANDARD_REQUIRED ON)

if(ENABLE_PROFILE MATCHES ON) 
  target_link_libraries(KaHyPar ${PROFILE_FLAGS})
  target_link_libraries(KaHyParShm ${PROFILE_FLAGS})
endif()

