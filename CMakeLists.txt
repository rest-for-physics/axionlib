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

set( external_libs ${external_libs} -lmpfr )

COMPILELIB("")

INSTALL(DIRECTORY ./data/
    DESTINATION ./data/axion/
    COMPONENT install
    )
