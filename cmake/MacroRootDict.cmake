if (APPLE)
    set(LD_LIBRARY_PATH_VAR DYLD_LIBRARY_PATH)
else ()
    set(LD_LIBRARY_PATH_VAR LD_LIBRARY_PATH)
endif ()
set(LD_LIBRARY_PATH_CONTENTS $ENV{${LD_LIBRARY_PATH_VAR}})
# MESSAGE( STATUS "LD_LIBRARY_PATH_CONTENTS: ${LD_LIBRARY_PATH_CONTENTS}" )

set(ROOT_CINT_WRAPPER
    ${LD_LIBRARY_PATH_VAR}=${ROOT_LIBRARY_DIR}:${LD_LIBRARY_PATH_CONTENTS}
    ${ROOTCINT_EXECUTABLE})

if (NOT DEFINED ROOT_DICT_OUTPUT_DIR)
    set(ROOT_DICT_OUTPUT_DIR "${PROJECT_BINARY_DIR}/rootdict")
endif ()

# clean generated header files with 'make clean'
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
                                    "${ROOT_DICT_OUTPUT_DIR}")

if (NOT ROOT_FIND_QUIETLY)
    message(
        STATUS "Check for ROOT_DICT_OUTPUT_DIR: ${PROJECT_BINARY_DIR}/rootdict")
    message(
        STATUS
            "Check for ROOT_DICT_CINT_DEFINITIONS: ${ROOT_DICT_CINT_DEFINITIONS}"
    )
endif ()

# ============================================================================
# helper macro to prepare input headers for GEN_ROOT_DICT_SOURCES sorts
# LinkDef.h to be the last header (required by rootcint)
#
# arguments: input_dir - directory to search for headers matching *.hh
#
# returns: ROOT_DICT_INPUT_HEADERS - all header files found in input_dir with
# ${input_dir}_LinkDef.h as the last header (if found)
#
# ----------------------------------------------------------------------------
macro (PREPARE_ROOT_DICT_HEADERS _input_dir)

    file(GLOB ROOT_DICT_INPUT_HEADERS "${_input_dir}/*.h")
    file(GLOB _linkdef_hdr "${_input_dir}/LinkDef.h")

    # LIST( FIND ROOT_DICT_INPUT_HEADERS ${_linkdef_hdr} _aux ) IF( ${_aux}
    # EQUAL 0 OR ${_aux} GREATER 0 ) LIST( REMOVE_ITEM ROOT_DICT_INPUT_HEADERS
    # "${_linkdef_hdr}" ) LIST( APPEND ROOT_DICT_INPUT_HEADERS "${_linkdef_hdr}"
    # ) ENDIF()

    if (_linkdef_hdr)
        list(REMOVE_ITEM ROOT_DICT_INPUT_HEADERS "${_linkdef_hdr}")
        list(APPEND ROOT_DICT_INPUT_HEADERS "${_linkdef_hdr}")
    endif ()

    # MESSAGE( STATUS "ROOT_DICT_INPUT_HEADERS: ${ROOT_DICT_INPUT_HEADERS}" )

endmacro (PREPARE_ROOT_DICT_HEADERS)

# ============================================================================
# helper macro to generate Linkdef.h files for rootcint
#
# arguments: namespace - prefix used for creating header <namespace>_Linkdef.h
# ARGN      - list of sources to be used for generating Linkdef.h
#
# returns: ROOT_DICT_INPUT_HEADERS - all header files + <namespace>_LinkDef.h in
# the correct order to be used by macro GEN_ROOT_DICT_SOURCES
#
# ----------------------------------------------------------------------------
macro (GEN_ROOT_DICT_LINKDEF_HEADER _namespace)

    set(_input_headers ${ARGN})
    set(_linkdef_header "${ROOT_DICT_OUTPUT_DIR}/${_namespace}_Linkdef.h")

    foreach (_header ${_input_headers})
        set(${_namespace}_file_contents
            "${${_namespace}_file_contents}\\#ifdef __CINT__\\\\n")
        set(${_namespace}_file_contents
            "${${_namespace}_file_contents}\\#pragma link off all globals\\;\\\\n"
        )
        set(${_namespace}_file_contents
            "${${_namespace}_file_contents}\\#pragma link off all classes\\;\\\\n"
        )
        set(${_namespace}_file_contents
            "${${_namespace}_file_contents}\\#pragma link off all functions\\;\\\\n"
        )
        set(${_namespace}_file_contents
            "${${_namespace}_file_contents}\\#pragma link C++ nestedclasses\\;\\\\n"
        )
        set(${_namespace}_file_contents
            "${${_namespace}_file_contents}\\#pragma link C++ nestedclasses\\;\\\\n"
        )
        set(${_namespace}_file_contents
            "${${_namespace}_file_contents}\\#pragma link C++ class ${_namespace}\\+\\;\\\\n"
        )
        set(${_namespace}_file_contents
            "${${_namespace}_file_contents}\\#endif\\\\n")
    endforeach ()

    add_custom_command(
        OUTPUT ${_linkdef_header}
        COMMAND mkdir -p ${ROOT_DICT_OUTPUT_DIR}
        COMMAND printf "${${_namespace}_file_contents}" > ${_linkdef_header}
        DEPENDS ${_input_headers}
        COMMENT "generating: ${_linkdef_header}")

    set(ROOT_DICT_INPUT_HEADERS ${_input_headers} ${_linkdef_header})

endmacro ()

# ============================================================================
# macro for generating root dict sources with rootcint
#
# arguments: dict_src_filename - filename of the dictionary source (to be
# generated)
#
# requires following variables: ROOT_DICT_INPUT_HEADERS - list of headers needed
# to generate dict source * if $LinkDef.h is in the list it must be at the end
# !! ROOT_DICT_INCLUDE_DIRS - list of include dirs to pass to rootcint -I..
# ROOT_DICT_CINT_DEFINITIONS - extra definitions to pass to rootcint
# ROOT_DICT_OUTPUT_DIR - where dictionary source should be generated
#
# returns: ROOT_DICT_OUTPUT_SOURCES - list containing generated source and other
# previously generated sources

# ----------------------------------------------------------------------------
macro (GEN_ROOT_DICT_SOURCE _dict_src_filename)

    set(_input_depend ${ARGN})
    # TODO check for ROOT_CINT_EXECUTABLE

    # need to prefix all include dirs with -I
    set(_dict_includes)
    foreach (_inc ${ROOT_DICT_INCLUDE_DIRS})
        set(_dict_includes "${_dict_includes}\t-I${_inc}") # fg: the \t fixes a
                                                           # wired string
                                                           # expansion
        # SET( _dict_includes ${_dict_includes} -I${_inc} )
    endforeach ()

    # We modify the list of headers to be given to ROOTCINT command. We must
    # remove/clean the full path from the main header
    list(GET ROOT_DICT_INPUT_HEADERS 0 MAIN_HEADER)
    get_filename_component(MAIN_HEADER_CLEAN ${MAIN_HEADER} NAME)
    list(GET ROOT_DICT_INPUT_HEADERS 1 LINKDEF_HEADER)
    set(ROOT_DICT_INPUT_HEADERS_CLEAN ${MAIN_HEADER_CLEAN} ${LINKDEF_HEADER})

    string(REPLACE "/" "_" _dict_src_filename_nosc ${_dict_src_filename})
    set(_dict_src_file ${ROOT_DICT_OUTPUT_DIR}/${_dict_src_filename_nosc})
    string(REGEX REPLACE "^(.*)\\.(.*)$" "\\1.h" _dict_hdr_file
                         "${_dict_src_file}")
    add_custom_command(
        OUTPUT ${_dict_src_file}
        COMMAND mkdir -p ${ROOT_DICT_OUTPUT_DIR}
        COMMAND ${ROOT_CINT_WRAPPER} -f "${_dict_src_file}" ${_dict_includes}
                ${ROOT_DICT_INPUT_HEADERS_CLEAN}
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        DEPENDS ${ROOT_DICT_INPUT_HEADERS} ${_input_depend}
        COMMENT "generating: ${_dict_src_file} with ${ROOT_DICT_INPUT_HEADERS}")
    list(APPEND ROOT_DICT_OUTPUT_SOURCES ${_dict_src_file})

endmacro ()

# for backwards compatibility
macro (GEN_ROOT_DICT_SOURCES _dict_src_filename)
    # MESSAGE( "USING DEPRECATED GEN_ROOT_DICT_SOURCES. PLEASE USE
    # GEN_ROOT_DICT_SOURCE instead." )
    set(ROOT_DICT_OUTPUT_SOURCES)
    gen_root_dict_source(${_dict_src_filename})
endmacro ()
# ============================================================================
