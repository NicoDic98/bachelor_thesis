set(DOXYGEN_IN "${PROJECT_SOURCE_DIR}/${DOCUMENTATION_DIR_NAME}/Doxyfile.in")
set(DOXYGEN_OUT "${PROJECT_BINARY_DIR}/${DOCUMENTATION_DIR_NAME}/Doxyfile")
configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

message("DOXYGEN_EXECUTABLE=${DOXYGEN_EXECUTABLE}")
add_custom_target(docs_doxygen ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        COMMENT "Generating documentation using Doxygen"
        VERBATIM)