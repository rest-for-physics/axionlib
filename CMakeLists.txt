#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
include(ExternalProject)
set(CMAKE_MACOSX_RPATH 1)

set( LibraryVersion "1.1" )
add_definitions(-DLIBRARY_VERSION="${LibraryVersion}")

find_path( MPFR_FOUND mpfr.h )
if( MPFR_FOUND MATCHES "MPFR_FOUND-NOTFOUND" )
	if( EXISTS ${MPFR_PATH}/include/mpfr.h )
		set( MPFR_FOUND FOUND )
	endif()
endif()

if( MPFR_FOUND MATCHES "MPFR_FOUND-NOTFOUND" )
	message( FATAL_ERROR "MPFR library was not found. Adding the location to the installation path may solve this problem.\n Use -DMPFR_PATH=/path/to/mpfr/" )
else()
	set( external_include_dirs ${external_include_dirs} ${MPFR_PATH}/include )
	link_directories( ${MPFR_PATH}/lib )
endif()

#-------------------------------------------------------------------------------------------------------
# Add the SolarAxionFlux library submodule as an external project.

if (${REST_SOLAXFLUX} MATCHES "ON")
    ## Integration of REST_SolaxFlux comes now as an option.
    ## It causes problems when active development happens at that library
    ## It should be fixed a versioning system. Probably through choosing the
    ## right submodule commit that works with RestAxionLib.
    ##
    ## For the moment I will just encapsulate SolaxFlux library to be
    ## able to compile and use a dummy/standard solar axion flux generation.
    add_definitions(-DUSE_SolaxFlux)

    set(SOLAXFLUX_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/solarAxionFlux")
    message( STATUS "SOLAXFLUX_DIR: ${SOLAXFLUX_DIR}" )

    ExternalProject_Add(solaxflux
                        SOURCE_DIR ${SOLAXFLUX_DIR}
                        GIT_SUBMODULES ""
                        GIT_CONFIG advice.detachedHead=False
                        CMAKE_CACHE_ARGS "-DPYTHON_SUPPORT:STRING=OFF;-DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER}")

    link_directories( ${SOLAXFLUX_DIR}/lib )
    set( external_include_dirs ${external_include_dirs} ${SOLAXFLUX_DIR}/include )
    set( external_libs ${external_libs} -laxionflux)
endif (${REST_SOLAXFLUX} MATCHES "ON")

set( external_libs ${external_libs} -lmpfr )

COMPILELIB("")

INSTALL(DIRECTORY ./data/
    DESTINATION ./data/axion/
    COMPONENT install
    )

INSTALL(DIRECTORY ./examples/
    DESTINATION ./examples/axion/
    COMPONENT install
    )

file(GLOB_RECURSE MAC "${CMAKE_CURRENT_SOURCE_DIR}/macros/*")
INSTALL(FILES ${MAC} DESTINATION ./macros/axion)
