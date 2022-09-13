Installation for Developers
===========================

In order to start developing, clone the QMAP repository using

    .. code-block:: console

        $ git clone --recurse-submodules -j8 https://github.com/cda-tum/qmap

Note the :code:`--recurse-submodules` flag. It is required to also clone all the required submodules. If you happen to forget passing the flag on your initial clone, you can initialize all the submodules by executing :code:`git submodule update --init --recursive` in the main project directory.

A C++ compiler supporting *C++17* and a minimum CMake version of *3.14* is required to build the project.

:code:`boost/program_options >= 1.50` is required for building the commandline applications of the mapping tool.

In order to build the exact mapping tool and for the Python bindings to work, the SMT Solver `Z3 >= 4.8.15 <https://github.com/Z3Prover/z3>`_  has to be installed and the dynamic linker has to be able to find the library. This can be accomplished in a multitude of ways:

- Under Ubuntu 20.04 and newer: :code:`sudo apt-get install libz3-dev`
- Under macOS: :code:`brew install z3`
- Alternatively: :code:`pip install z3-solver` and then append the corresponding path to the library path (:code:`LD_LIBRARY_PATH` under Linux, :code:`DYLD_LIBRARY_PATH` under macOS), e.g. via

    .. code-block:: console

        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(python -c "import z3; print(z3.__path__[0]+'/lib')")

- Download pre-built binaries from https://github.com/Z3Prover/z3/releases and copy the files to the respective system directories
- Build Z3 from source and install it to the system

.. note::
    We noticed some issues when compiling with Microsoft's *MSCV* compiler toolchain. If you want to start development on this project under Windows, consider using the *clang* compiler toolchain. A detailed description of how to set this up can be found `here <https://docs.microsoft.com/en-us/cpp/build/clang-support-msbuild?view=msvc-160>`_.

Working on the core C++ library
###############################

Our projects use *CMake* as the main build configuration tool.
Building a project using CMake is a two-stage process.
First, CMake needs to be *configured* by calling

    .. code-block:: console

        $ cmake -S . -B build -DCMAKE_BUILD_TYPE=Release

This tells CMake to search the current directory :code:`.` (passed via :code:`-S`) for a :code:`CMakeLists.txt` file and process it into a directory :code:`build` (passed via :code:`-B`).
The flag :code:`-DCMAKE_BUILD_TYPE=Release` tells CMake to configure a *Release* build (as opposed to, e.g., a *Debug* build).

After configuring with CMake, the project can be built by calling

    .. code-block:: console

        $ cmake --build build --config Release

This tries to build the project in the :code:`build` directory (passed via :code:`--build`).
Some operating systems and developer environments explicitly require a configuration to be set, which is why the :code:`--config` flag is also passed to the build command. The flag :code:`--parallel <NUMBER_OF_THREADS>` may be added to trigger a parallel build.

Building the project this way generates

- the heuristic library :code:`libqmap_heuristic_lib.a` (Unix) / :code:`qmap_heuristic_lib.lib` (Windows) in the :code:`build/src` directory
- the heuristic mapper commandline executable :code:`qmap_heuristic` in the :code:`build/apps` directory (only available if Boost is found)
- a test executable :code:`qmap_heuristic_test` containing a small set of unit tests for the heuristic mapper in the :code:`build/test` directory
- the exact library :code:`libqmap_exact_lib.a` (Unix) / :code:`qmap_exact_lib.lib` (Windows) in the :code:`build/src` directory (only available if Z3 is found)
- the exact mapper commandline executable :code:`qmap_exact` in the :code:`build/apps` directory (only available if Boost and Z3 is found)
- a test executable :code:`qmap_exact_test` containing a small set of unit tests for the exact mapper in the :code:`build/test` directory (only available if Z3 is found)

Working on the Python module
############################

The :code:`mqt.qmap` Python module can be conveniently built locally by calling

    .. code-block:: console

        (venv) $ pip install --editable .[dev]

The :code:`--editable` flag ensures that changes in the Python code are instantly available without re-running the command. The :code:`[dev]` extra makes sure that all dependencies for running the Python tests and building the documentation are available.

.. note::
    When using the :code:`zsh` shell it might be necessary to add double quotes around the :code:`.[dev]` part of the command.

`Pybind11 <https://pybind11.readthedocs.io/>`_ is used for providing bindings of the C++ core library to Python (see `bindings.cpp <https://github.com/cda-tum/qmap/tree/main/mqt/qmap/bindings.cpp>`_).
If parts of the C++ code have been changed, the above command has to be run again to make the changes visible in Python.

Running tests
#############

C++ core library
----------------

The C++ part of the code base is tested by unit tests using the `googletest <https://google.github.io/googletest/primer.html>`_ framework.
The corresponding test files can be found in the :code:`test` directory. In order to build the tests, CMake first has to be configured with the :code:`-DBUILD_QMAP_TESTS=ON` flag, i.e.,

    .. code-block:: console

        $ cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_QMAP_TESTS=ON

Then, the test executables :code:`qmap_heuristic_test` and :code:`qmap_exact_test` is built in the :code:`build/test` directory by calling

    .. code-block:: console

        $ cmake --build build --config Release --target qmap_heuristic_test
        $ cmake --build build --config Release --target qmap_exact_test

From there, the tests can be started by simply calling

    .. code-block:: console

        [.../build/test] $ ./qmap_heuristic_test
        [.../build/test] $ ./qmap_exact_test

.. note::
    The test executable :code:`qmap_exact_test` will only be built if the Z3 solver is available during building.

Python interface and functionality
----------------------------------

The Python part of the code base is tested by unit tests using the `pytest <https://docs.pytest.org/en/latest/>`_ framework.
The corresponding test files can be found in the :code:`test/python` directory.
To start the tests, simply call

    .. code-block:: console

        (venv) $ python -m pytest ./test/python

.. note::
    If you haven't already installed the package with the :code:`[dev]` extra as demonstrated above, the necessary dependencies for running the Python tests can be installed by calling

        .. code-block:: console

            (venv) $ pip install --editable .[test]

Building the documentation
##########################

Building this documentation locally is as easy as calling

    .. code-block:: console

        (venv) [.../docs] $ make clean && make html

The resulting HTML documentation (:code:`index.html`) can be found in the :code:`docs/build/html` directory.

.. note::
    If you haven't already installed the package with the :code:`[dev]` extra as demonstrated above, the necessary dependencies for building the documentation can be installed by calling

        .. code-block:: console

            (venv) $ pip install --editable .[docs]
