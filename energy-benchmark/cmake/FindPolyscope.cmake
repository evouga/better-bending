# - Try to find the Polyscope library
# Once done this will define
#
#  POLYSCOPE_FOUND - system has POLYSCOPE 
#  POLYSCOPE_INCLUDE_DIR - **the** POLYSCOPE include directory
if(TARGET polyscope)
	if(NOT POLYSCOPE_INCLUDE_DIR)
		find_path(POLYSCOPE_INCLUDE_DIR polyscope/polyscope.h
			HINTS
				ENV POLYSCOPE_DIR
			PATHS
				${CMAKE_SOURCE_DIR}/../..
				${CMAKE_SOURCE_DIR}/..
				${CMAKE_SOURCE_DIR}
				${CMAKE_SOURCE_DIR}/polyscope
				${CMAKE_SOURCE_DIR}/../polyscope
				${CMAKE_SOURCE_DIR}/../tools/polyscope
				${CMAKE_SOURCE_DIR}/../../polyscope
				/usr
				/usr/local
			PATH_SUFFIXES include
		)
	endif()
	set(POLYSCOPE_FOUND TRUE)
	set(Polyscope_FOUND TRUE)
	return()
endif()

if(POLYSCOPE_FOUND)
	return()
endif()

find_path(POLYSCOPE_INCLUDE_DIR polyscope/polyscope.h
    HINTS
        ENV POLYSCOPE_DIR
    PATHS
        ${CMAKE_SOURCE_DIR}/../..
        ${CMAKE_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/polyscope
        ${CMAKE_SOURCE_DIR}/../polyscope
        ${CMAKE_SOURCE_DIR}/../tools/polyscope
        ${CMAKE_SOURCE_DIR}/../../polyscope
        /usr
        /usr/local
    PATH_SUFFIXES include
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Polyscope 
    "\npolyscope not found"
    POLYSCOPE_INCLUDE_DIR)
mark_as_advanced(POLYSCOPE_INCLUDE_DIR)

if(NOT TARGET polyscope)
	add_subdirectory("${POLYSCOPE_INCLUDE_DIR}/../" "polyscope")
endif()