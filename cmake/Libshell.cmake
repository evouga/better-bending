# Build the libshell static library from the submodule without running
# libshell/CMakeLists.txt. The upstream libshell build always configures its
# example target, which FetchContent's a second copy of polyscope and conflicts
# with the polyscope submodule used by the experiment projects.
function(better_bending_add_libshell)
	if(TARGET libshell)
		return()
	endif()

	list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/libshell/cmake")
	include(libigl)

	file(GLOB LIBSHELL_SOURCES
		"${CMAKE_SOURCE_DIR}/libshell/src/*.cpp"
		"${CMAKE_SOURCE_DIR}/libshell/src/SecondFundamentalForm/*.cpp"
		"${CMAKE_SOURCE_DIR}/libshell/src/MaterialModel/*.cpp"
	)
	add_library(libshell STATIC ${LIBSHELL_SOURCES})
	target_include_directories(libshell PUBLIC "${CMAKE_SOURCE_DIR}/libshell/include")
	target_link_libraries(libshell PUBLIC Eigen3::Eigen)
	set_target_properties(libshell PROPERTIES CXX_STANDARD 11)

	if(MSVC)
		target_compile_definitions(libshell PRIVATE EIGEN_STRONG_INLINE=inline)
		target_compile_options(libshell PRIVATE /bigobj)
	endif()
endfunction()
