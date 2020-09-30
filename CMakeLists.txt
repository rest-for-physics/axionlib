set( THIS_LIBRARY "RestAxion" )
set( LibraryVersion "1.0" )
add_definitions(-DLIBRARY_VERSION="${LibraryVersion}")
set( external_libs ${external_libs} -lmpfr )

set(contents)

COMPILEDIR(RestAxion)
