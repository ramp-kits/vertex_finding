cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
cmake_policy(SET CMP0042 NEW)

project(baselineAPP C CXX)


set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS 1)

set(CMAKE_BUILD_TYPE RelWithDebInfo)


list(APPEND CMAKE_PREFIX_PATH ${CMAKE_INSTALL_PREFIX}) # for find_package
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -march=native")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG")

set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 



find_package(Boost)
add_subdirectory(baseline)
add_subdirectory(PatPV)
include_directories(PatPV/include ${Boost_INCLUDE_DIRS})
include_directories(baseline/baseline)
file(GLOB ramp_sources "${CMAKE_CURRENT_SOURCE_DIR}/baseline/baseline.cpp")
add_executable(ramp_test ${ramp_sources} )
target_link_libraries(ramp_test  PUBLIC RAMP x86PatPV  ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

