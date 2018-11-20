include_directories(../PatPV/include)

find_package(Boost)
include_directories(${Boost_INCLUDE_DIRS})

add_library(marianRAMP  SHARED
  readData.cpp
  baseline.cpp
)
target_link_libraries( marianRAMP PUBLIC marianPatPV  ${ROOT_LIBRARIES} ${BOOST_LIBRARIES} )
