# Build Python bindings
include(FetchContent)
FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11
        GIT_TAG        v2.2.4
)

FetchContent_GetProperties(pybind11)
if(NOT pybind11_POPULATED)
    FetchContent_Populate(pybind11)
    add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
endif()
set(PYBIND11_CPP_STANDARD -std=c++17)

include_directories(pybind11/include)

set(STUKAPY_SOURCES
        src/stukapy.cpp
        src/lp.h
        src/lp.cpp
        src/qp.h
        src/qp.cpp
        src/dlp.h
        src/dlp.cpp
        src/stochastic.h
        src/stochastic.cpp
        src/callbacks.h
        src/callbacks.cpp)

pybind11_add_module(stukapy ${STUKAPY_SOURCES})
target_link_libraries(stukapy PRIVATE stuka)
set_target_properties(stukapy PROPERTIES CXX_VISIBILITY_PRESET hidden)