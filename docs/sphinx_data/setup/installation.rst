===================================
Installation
===================================

This library takes advantage of two possible build systems - Bazel_ (recommended) or CMake_. It additional has the following dependencies.

- C++17 compiler (ex. gcc >= 7 or clang >= 5)
- OpenBLAS_ (required for Bazel_ and optional for CMake_ if another LAPACK and BLAS implementation is available)
- SuiteSparse_
- HDF5_
- Gurobi_ (optional)

Downloading Stuka
===================================

Open a terminal and download Stuka with

.. code-block:: bash

    git clone https://github.com/amadavan/Stuka.git

Building with Bazel
===================================

The Bazel_ build system offers significant advantages over the CMake_ system. In particular, it creates hermetic builds. That is, the constructed files are independent of any other sources. Additionally, the build system will automatically construct the dependencies.

.. code-block:: bash

    bazel build //:stuka

Optional arguments include :code:`--config=gurobi` if you have gurobi installed or :code:`--config=mkl` if you wish to use MKL's implementation of BLAS and LAPACK. Note that this requires Intel's MKL to be located and available on the path. If you a local OpenBLAS installation, or HDF5 installed locally then you may not want Bazel_ to build all of these. Stuka can use local libraries, by setting the environment variable :code:`STUKA_SYSTEM_LIBS={SYSTEM_LIBRARIES}`, where :code:`{SYSTEM_LIBRARIES}` is a comma-delimited list of system libraries to use.

.. warning::
    The first build of Stuka will take a long time to build if an OpenBLAS system library is not provided. This is due to the fact that Stuka must build the entire OpenBLAS library from scratch. Unfortunately, the Bazel build cannot currently use any other BLAS or LAPACK implementations.

Building with CMake
===================================
The CMake_ build system requires manual installation of the dependencies.

Installing SuiteSparse
---------------------------------
The scripts folder of Stuka contains a way to build SuiteSparse_. The scripts takes advantage of the provided Makefile which only constructs the SuiteSparse_ targets that are required for the library. The following command will ensure that the script has executable rights and then install SuiteSparse_ in the lib folder.

.. code-block:: bash

    chmod 755 ./scripts/install_suitesparse.sh
    ./scripts/install_suitesparse.sh

Installing HDF5
---------------------------------
The HDF5 library can be installed directly from source by running the following commands on a UNIX-like system,

.. code-block:: bash

    wget HDF5_link
    ./configure
    make
    make install

Building Stuka
---------------------------------
Create a new directory called :file:`build` in which Stuka will be built. Then run the following commands,

.. code-block:: bash

    cmake .. -DCMAKE_BUILD_TYPE=Release
    make -j4 stuka

.. _OpenBLAS: https://www.openblas.net/
.. _SuiteSparse: http://faculty.cse.tamu.edu/davis/suitesparse.html
.. _HDF5: https://www.hdfgroup.org/solutions/hdf5/
.. _Gurobi: http://www.gurobi.com/downloads/download-center
.. _Bazel: https://bazel.build/
.. _CMake: https://cmake.org/