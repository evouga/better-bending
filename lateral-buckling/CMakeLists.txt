cmake_minimum_required(VERSION 4.0)
cmake_policy(SET CMP0169 OLD)
project(lateral-buckling)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()


include(libigl)
igl_include(restricted triangle)

find_package(Polyscope REQUIRED)
find_package(Libshell REQUIRED)


file(GLOB SRCFILES src/*.cpp src/*.h src/*.hpp)	  
file(GLOB VIEWERFILES viewer/*.cpp viewer/*.h src/Visualization.h src/Visualization.cpp)

add_executable(${PROJECT_NAME} ${SRCFILES})
target_include_directories(${PROJECT_NAME} PUBLIC ${EIGEN3_INCLUDE_DIR} ${LIBSHELL_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} ${LIBSHELL_LIB_DIR} igl_restricted::triangle igl::core polyscope)

add_executable(Viewer ${VIEWERFILES})
target_include_directories(Viewer PUBLIC ${EIGEN3_INCLUDE_DIR} ${LIBSHELL_INCLUDE_DIR})
target_link_libraries(Viewer ${LIBSHELL_LIB_DIR} igl::core polyscope)

# Use C++14
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 14)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD_REQUIRED ON)
