include_directories(../PatPV/include)

find_package(Boost)
include_directories(${Boost_INCLUDE_DIRS})

add_library(RAMP  SHARED
  readData.cpp
  baseline.cpp
)
target_link_libraries( RAMP PUBLIC x86PatPV  ${ROOT_LIBRARIES} ${BOOST_LIBRARIES} )
