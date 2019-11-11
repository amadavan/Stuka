if (DOXYGEN_FOUND)

    set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
    set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    #Replace variables inside @@ with the current values
    configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

    add_custom_target(doxygen
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating docs")

    get_target_property(DOXYGEN_OUTPUT_DIR doxygen BINARY_DIR)

    if (SPHINX_FOUND)

        set(DOXYGEN_OUTPUT_DIR "${DOXYGEN_OUTPUT_DIR}/xml")

        add_custom_target(sphinx ALL COMMAND
                ${SPHINX_EXECUTABLE} -b html
                -Dbreathe_projects.stuka=${DOXYGEN_OUTPUT_DIR}
                ${CMAKE_CURRENT_SOURCE_DIR}/sphinx_data ${CMAKE_CURRENT_BINARY_DIR}
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                COMMENT "Generating documentation with Sphinx")
        add_dependencies(sphinx doxygen)

    endif (SPHINX_FOUND)

endif (DOXYGEN_FOUND)
