#----------------------------------------------------------------------------
# Install examples to the 'examples' directory and link this library
#----------------------------------------------------------------------------
add_executable(run_first_example src/FirstExample.cxx)
target_link_libraries(run_first_example LINK_PUBLIC ${THIS_LIBRARY})
set_target_properties(run_first_example
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/examples/"
)
