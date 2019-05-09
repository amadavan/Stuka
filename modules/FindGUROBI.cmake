#### Taken from http://www.openflipper.org/svnrepo/CoMISo/trunk/CoMISo/cmake/FindGUROBI.cmake


# - Try to find GUROBI
# Once done this will define
#  GUROBI_FOUND - System has Gurobi
#  GUROBI_INCLUDE_DIRS - The Gurobi include directories
#  GUROBI_LIBRARIES - The libraries needed to use Gurobi

if (GUROBI_INCLUDE_DIR)
    # in cache already
    set(GUROBI_FOUND TRUE)
    set(GUROBI_INCLUDE_DIRS "${GUROBI_INCLUDE_DIR}")
    set(GUROBI_LIBRARIES "${GUROBI_LIBRARY};${GUROBI_CXX_LIBRARY}")
else (GUROBI_INCLUDE_DIR)

    find_program(LSB_RELEASE_EXEC lsb_release)
    execute_process(COMMAND ${LSB_RELEASE_EXEC} -is
            OUTPUT_VARIABLE LSB_RELEASE_ID_SHORT
            OUTPUT_STRIP_TRAILING_WHITESPACE
            )
    if (${LSB_RELEASE_ID} MATCHES "Ubuntu")
        set(GUROBI_CXX_LIBNAME gurobi_c++5.2 gurobi_c++4.2)
    endif ()

    find_path(GUROBI_INCLUDE_DIR
            NAMES gurobi_c++.h
            PATHS "$ENV{GUROBI_HOME}/include"
            "/Library/gurobi811/mac64/include"
            "C:\\libs\\gurobi811\\include"
            "/Library/gurobi800/mac64/include"
            "C:\\libs\\gurobi800\\include"
            "/Library/gurobi752/mac64/include"
            "C:\\libs\\gurobi752\\include"
            )

    find_library(GUROBI_LIBRARY
            NAMES gurobi
            gurobi81
            gurobi80
            gurobi75
            PATHS "$ENV{GUROBI_HOME}/lib"
            "/Library/gurobi811/mac64/lib"
            "C:\\libs\\gurobi811\\lib"
            "/Library/gurobi800/mac64/lib"
            "C:\\libs\\gurobi800\\lib"
            "/Library/gurobi752/mac64/lib"
            "C:\\libs\\gurobi752\\lib"
            )

    find_library(GUROBI_CXX_LIBRARY
            NAMES ${GUROBI_CXX_LIBNAME} gurobi_c++
            PATHS "$ENV{GUROBI_HOME}/lib"
            "/Library/gurobi811/mac64/lib"
            "C:\\libs\\gurobi811\\lib"
            "/Library/gurobi800/mac64/lib"
            "C:\\libs\\gurobi800\\lib"
            "/Library/gurobi752/mac64/lib"
            "C:\\libs\\gurobi752\\lib"
            )

    set(GUROBI_INCLUDE_DIRS "${GUROBI_INCLUDE_DIR}")
    set(GUROBI_LIBRARIES "${GUROBI_CXX_LIBRARY};${GUROBI_LIBRARY}")

    # use c++ headers as default
    # set(GUROBI_COMPILER_FLAGS "-DIL_STD" CACHE STRING "Gurobi Compiler Flags")

    include(FindPackageHandleStandardArgs)
    # handle the QUIETLY and REQUIRED arguments and set LIBCPLEX_FOUND to TRUE
    # if all listed variables are TRUE
    find_package_handle_standard_args(GUROBI DEFAULT_MSG
            GUROBI_LIBRARY GUROBI_CXX_LIBRARY GUROBI_INCLUDE_DIR)

    mark_as_advanced(GUROBI_INCLUDE_DIR GUROBI_LIBRARY GUROBI_CXX_LIBRARY)

endif (GUROBI_INCLUDE_DIR)