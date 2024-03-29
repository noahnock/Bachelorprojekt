cmake_minimum_required(VERSION 3.19)
project(Bachelorprojekt)

set(CMAKE_CXX_STANDARD 20)

add_library(rtree
        Rtree/Rtree.h
        Rtree/Rtree.cpp
        Rtree/RtreeBuild.cpp
        Rtree/RtreeFileReader.cpp
        Rtree/RtreeFileReader.h
        ExternalSorting.cpp
        Rtree/RtreeBasicGeometry.h
        Rtree/RtreeNode.cpp
        Rtree/RtreeNode.h
        Rtree/RtreeSorter.cpp
        Rtree/RtreeSorter.h)
target_link_libraries(rtree Boost::serialization Boost::program_options Boost::iostreams)

add_executable(Bachelorprojekt main.cpp)
target_link_libraries(Bachelorprojekt ${Boost_LIBRARIES} Boost::serialization Boost::program_options Boost::iostreams rtree)


find_package(Boost 1.76 COMPONENTS iostreams serialization program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Wall -Wextra")
# Enable the specification of additional warnings manually from the commandline
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ADDITIONAL_COMPILER_FLAGS}")


if (${PERFTOOLS_PROFILER})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lprofiler")
    message(STATUS "Adding -lprofiler (make sure your have google-perftools installed.)")
endif ()

if (${ALLOW_SHUTDOWN})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DALLOW_SHUTDOWN")
    message(STATUS "Adding -DALLOW_SHUTDOWN")
endif ()


set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")


set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}")