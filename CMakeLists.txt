set( LibraryVersion "1.0" )
add_definitions(-DLIBRARY_VERSION="${LibraryVersion}")

set( external_libs ${external_libs} -lmpfr )

COMPILELIB("")

INSTALL(DIRECTORY ./data/
    DESTINATION ./data/axion/
    COMPONENT install
    )
