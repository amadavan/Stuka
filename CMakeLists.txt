cmake_minimum_required(VERSION 3.0)

project(Stuka)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS_DEBUG "-g -fsanitize=address")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib CACHE PATH "Shared-object library path")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib CACHE PATH "Static library path")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin CACHE PATH "Runtime path")

set(CMAKE_INCLUDE_PATH ${CMAKE_SOURCE_DIR}/include)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/modules")
list(APPEND INCLUDE_INSTALL_DIR "${CMAKE_SOURCE_DIR}/include")
list(APPEND LIB_INSTALL_DIR "${CMAKE_SOURCE_DIR}/lib")

# Required packages
find_package(EIGEN REQUIRED)
find_package(CHOLMOD REQUIRED)
find_package(SPQR REQUIRED)

# Optional packages
find_package(GUROBI)
find_package(HDF5 COMPONENTS C)
find_package(DOXYGEN)
find_package(SPHINX)

add_subdirectory(src)
add_subdirectory(bindings)
add_subdirectory(docs)