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
    if (${LSB_RELEASE_ID_SHORT} MATCHES "Ubuntu")
        set(GUROBI_CXX_LIBNAME gurobi_g++5.2)
    else()
        set(GUROBI_CXX_LIBNAME gurobi_c++)
    endif ()

    # Get version and set the gurobi root directory
    if (DEFINED ENV{GUROBI_HOME})
        string(REGEX MATCH "[0-9]+" GUROBI_VERSION "$ENV{GUROBI_HOME}")
        set(GUROBI_DIR ${GUROBI_HOME})
    else()
        if(APPLE)
            file(GLOB dirs /Library/gurobi*)
            foreach(d in ${dirs})
                string(REGEX MATCH "[0-9]+" GUROBI_VERSION "${d}")
            endforeach(d)
            set(GUROBI_DIR /Library/gurobi${GUROBI_VERSION}/mac64)
        elseif(UNIX)
            file(GLOB dirs /opt/gurobi/gurobi*)
            foreach(d in ${dirs})
                string(REGEX MATCH "[0-9]+" GUROBI_VERSION "${d}")
            endforeach(d)
            set(GUROBI_DIR /opt/gurobi/gurobi${GUROBI_VERSION}/linux64)
        elseif(WIN32 OR MSVC)
            file(GLOB dirs "C:\\libs\\gurobi*")
            foreach(d in ${dirs})
                string(REGEX MATCH "[0-9]+" GUROBI_VERSION "${d}")
            endforeach()
            set(GUROBI_DIR "C:\\libs\\gurobi${GUROBI_VERSION}")
        endif()
    endif()
    string(SUBSTRING "${GUROBI_VERSION}" 0 2 GUROBI_VERSION_SHORT)

    message("Using GUROBI version ${GUROBI_VERSION} from ${GUROBI_DIR}")

    find_path(GUROBI_INCLUDE_DIR
            NAMES gurobi_c++.h
            PATHS "$ENV{GUROBI_HOME}/include"
            "${GUROBI_DIR}/include"
            "${GUROBI_DIR}\\include"
            )

    find_library(GUROBI_LIBRARY
            NAMES gurobi${GUROBI_VERSION_SHORT}
            PATHS "$ENV{GUROBI_HOME}/lib"
            "${GUROBI_DIR}/lib"
            "${GUROBI_DIR}\\lib"
            )

    find_library(GUROBI_CXX_LIBRARY
            NAMES ${GUROBI_CXX_LIBNAME}
            PATHS "$ENV{GUROBI_HOME}/lib"
            "${GUROBI_DIR}/lib"
            "${GUROBI_DIR}\\lib"
            )

    set(GUROBI_INCLUDE_DIRS "${GUROBI_INCLUDE_DIR}")
    set(GUROBI_LIBRARIES "${GUROBI_CXX_LIBRARY};${GUROBI_LIBRARY}")
    message(${GUROBI_CXX_LIBRARY})

    # use c++ headers as default
    # set(GUROBI_COMPILER_FLAGS "-DIL_STD" CACHE STRING "Gurobi Compiler Flags")

    include(FindPackageHandleStandardArgs)
    # handle the QUIETLY and REQUIRED arguments and set LIBCPLEX_FOUND to TRUE
    # if all listed variables are TRUE
    find_package_handle_standard_args(GUROBI DEFAULT_MSG
            GUROBI_LIBRARY GUROBI_CXX_LIBRARY GUROBI_INCLUDE_DIR)

    mark_as_advanced(GUROBI_INCLUDE_DIR GUROBI_LIBRARY GUROBI_CXX_LIBRARY)

endif (GUROBI_INCLUDE_DIR)