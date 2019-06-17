![Logo](assets/stuka-banner.png)
--------------------------------
[![Build Status](https://travis-ci.org/amadavan/Stuka.svg?branch=master)](https://travis-ci.org/amadavan/Stuka)

# Requirements
- C++14 compiler (gcc >= 5 or clang >= 3.5)
- [CMake](https://cmake.org/download/) (version >= 3.11) or [Bazel](https://github.com/bazelbuild/bazel/releases)
- [SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html)
- [Gurobi optimizer](http://www.gurobi.com/downloads/download-center)

# Installation
#### SuiteSparse
Before installing SuiteSparse, make sure you have a `blas` and `lapack`
implementation installed. For convenience a suitesparse installer is
available in the `scripts/` directory. Make sure you are in the root
`Stuka` directory and run the following command.
```bash
./scripts/install_suitesparse
```
If the scripts fails due to an inability to find `blas` or `lapack` you
may need to manually specify the location of your implementation by
editing the `scripts/install_suitesparse.sh` file with the following
```bash
make library BLAS=${BLAS_LOCATION}
```

#### Stuka
Stuka is built using the CMake build system. To create the library
simply run the following commands from the root `Stuka` directory.
```bash
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . --target stuka
```
This should create a shared library in the `build/lib` directory.

#### Python bindings
Building the python library is as simple running
```bash
python setup.py install
```
This should work assuming there are no issues in the previous step.

# Examples
#### C++ example
In order to test the code, a test script is available to be built.
Assuming you have build the library as described in the Installation
section, run the following.
```bash
cd build
cmake --build . --target tester
```
This will create a `tester` executable in the `build/bin` directory.
Running the tester will run into an error, but this is expected. In
the final example, Bender's Decomposition is unable to determine the
optimal solution due to degeneracy of the problem.

#### Python example
This will be added at a later date.