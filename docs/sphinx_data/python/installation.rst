Installation
=================================

While we recommend Bazel for the C++ build, the current Bazel rules do not work for the python bindings. Instead, the python bindings rely on the CMake build system. In order to install the bindings, simply run the following from Stuka's root directory.

.. code-block:: bash

    python setup.py install

This will install the :code`stukapy` package which you can then use directly. It's as simple as that!