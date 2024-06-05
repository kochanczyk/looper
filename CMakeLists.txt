cmake_minimum_required(VERSION 3.5)

project(Looper CXX)


##### Compiler

set(CMAKE_CXX_STANDARD_REQUIRED 14)
set(CMAKE_CXX_EXTENSIONS OFF)
#
set(CMAKE_CXX_FLAGS         "-Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG   "-O1 -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3 -mtune=native")


##### Dependencies

### Dlib (mandatory)
#
set(INTERNAL_DLIB_VERSION "19.9")
set(Dlib_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/contrib/dlib-${INTERNAL_DLIB_VERSION}/include")
include_directories(SYSTEM ${Dlib_INCLUDE_DIR})
add_definitions(-DDLIB_NO_GUI_SUPPORT -UDLIB_PNG_SUPPORT -UDLIB_JPEG_SUPPORT)

### Boost (optional)
#
if (EXISTS ${Boost_INCLUDE_DIR})
  message(STATUS "Note: Using boost from ${Boost_INCLUDE_DIR}")
else()
  set(INTERNAL_BOOST_VERSION "1.64.0")
  set(Boost_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/contrib/boost_${INTERNAL_BOOST_VERSION}/include")
  set(Boost_ADDITIONAL_VERSIONS _INTERNAL_BOOST_VERSION)
  message(STATUS "Note: Using internal boost from ${Boost_INCLUDE_DIR}")
endif()
#
find_package(Boost)
if(Boost_FOUND)
  include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
  add_definitions(-DCTHULHU_USES_BOOST=1)
else()
  add_definitions(-DCTHULHU_USES_BOOST=0)
endif()


##### Components

set(SRC_LIST Looper.cpp
             Settings.cpp
             Settings.hpp
             Trajectory.hpp
             TrajectoryStoragePolicies.hpp
             Chemistry.hpp
             ThreadPool.hpp
             common.hpp
             "${Dlib_INCLUDE_DIR}/dlib/all/source.cpp")

set_source_files_properties("${Dlib_INCLUDE_DIR}/dlib/all/source.cpp" PROPERTIES 
                            COMPILE_FLAGS "-Wno-pedantic -Wno-extra -Wno-all")
add_executable(Looper ${SRC_LIST})
target_link_libraries(Looper pthread)

