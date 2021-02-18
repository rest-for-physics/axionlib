#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
include(ExternalProject)
set(CMAKE_MACOSX_RPATH 1)

## --------------------------------------------------------------------------
## ---------- ONLY need to change THIS_LIBRARY and LibraryVersion -----------
## ---------- to adapt this dummy CMakeLists to a new library name.  --------
## ---------- The version value should be assigned to any metadata ----------
## ---------- structure inside the Initialize method              -----------
## ----------
## ---------- THIS_LIBRARY should contain `Rest` or `REST` in the -----------
## ---------- name to be recognized by TRestTools::GetListOfRestLibraries. --
## ---------- You might use any string you like in the version field --------
## ----------
set( THIS_LIBRARY "RestAxion" )
set( LibraryVersion "1.0" )
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
#-------------------------------------------------------------------------------------------------------

set( external_libs ${external_libs} -lmpfr )

COMPILELIB("")

INSTALL(DIRECTORY ./data/
    DESTINATION ./data/axion/
    COMPONENT install
    )
