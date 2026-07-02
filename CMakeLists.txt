cmake_minimum_required(VERSION 4.0)
cmake_policy(SET CMP0169 OLD)
project(better-bending)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

# Build shared dependencies first.
add_subdirectory(libshell)

find_path(POLYSCOPE_INCLUDE_DIR polyscope/polyscope.h
	PATHS ${CMAKE_CURRENT_SOURCE_DIR}/polyscope
	PATH_SUFFIXES include
	NO_DEFAULT_PATH
)
if(NOT POLYSCOPE_INCLUDE_DIR)
	message(FATAL_ERROR
		"polyscope not found at ${CMAKE_CURRENT_SOURCE_DIR}/polyscope. "
		"Run: git submodule update --init --recursive")
endif()
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/polyscope ${CMAKE_BINARY_DIR}/polyscope)

# Child projects use FindLibshell.cmake / FindPolyscope.cmake, which look for
# prebuilt installs by default. Point them at the in-tree targets instead.
set(LIBSHELL_FOUND TRUE CACHE BOOL "libshell built from submodule" FORCE)
set(LIBSHELL_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/libshell/include" CACHE PATH "" FORCE)
set(LIBSHELL_LIB_DIR libshell CACHE STRING "" FORCE)

set(POLYSCOPE_FOUND TRUE CACHE BOOL "polyscope built from submodule" FORCE)
set(POLYSCOPE_INCLUDE_DIR ${POLYSCOPE_INCLUDE_DIR} CACHE PATH "" FORCE)

add_subdirectory(energy-benchmark)
add_subdirectory(lateral-buckling)
