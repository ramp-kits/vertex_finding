include_directories(include)
include_directories(../cuda/velo/common/include)
include_directories(../cuda/patPV/include ${Boost_INCLUDE_DIRS})

add_library(marianPatPV SHARED
  ${CMAKE_CURRENT_SOURCE_DIR}/src/PVSeedTool.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/AdaptivePV3DFitter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/run_PatPV_CPU.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/PVChecker.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/patPV_Definitions.cuh
)

target_link_libraries(marianPatPV PUBLIC ${Boost_LIBRARIES} ${ROOT_LIBRARIES})
# target_link_libraries(x86PatPV)
