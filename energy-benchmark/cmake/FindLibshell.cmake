# - Try to find the LIBSHELL library
# Once done this will define
#
#  LIBSHELL_FOUND - system has LIBSHELL 
#  LIBSHELL_INCLUDE_DIR - **the** LIBSHELL include directory
#  LIBCOLLISIONDETECTION_LIB_DIR - the LIBSHELL library directory
if(LIBSHELL_FOUND)
    return()
endif()

find_path(LIBSHELL_INCLUDE_DIR ElasticShell.h
    HINTS
        ENV LIBSHELL_DIR
    PATHS
        ${CMAKE_SOURCE_DIR}/../..
        ${CMAKE_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/libshell
        ${CMAKE_SOURCE_DIR}/../libshell
        ${CMAKE_SOURCE_DIR}/../../libshell
        /usr
        /usr/local        
    PATH_SUFFIXES include
)

find_library(LIBSHELL_LIB_DIR libshell
    HINTS
        ENV LIBSHELL_DIR
    PATHS
        ${CMAKE_SOURCE_DIR}/../..
        ${CMAKE_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/libshell
        ${CMAKE_SOURCE_DIR}/../libshell
        ${CMAKE_SOURCE_DIR}/../../libshell
        /usr
        /usr/local        
    PATH_SUFFIXES lib 
)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Libshell
    "\nlibshell not found"
    LIBSHELL_INCLUDE_DIR LIBSHELL_LIB_DIR)
mark_as_advanced(LIBSHELL_INCLUDE_DIR LIBSHELL_LIB_DIR)

