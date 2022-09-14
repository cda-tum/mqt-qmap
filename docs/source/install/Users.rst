Installation for Users
======================

QMAP is mainly developed as a C++ library that builds upon `our quantum functionality representation (QFR) <https://github.com/cda-tum/qfr>`_.
In order to make the tool as accessible as possible, it comes with an easy-to-use Python interface.

We encourage installing QMAP via pip (preferably in a `virtual environment <https://docs.python.org/3/library/venv.html>`_):

    .. code-block:: console

        (venv) $ pip install mqt.qmap

In most practical cases (under 64-bit Linux, MacOS incl. Apple Silicon, and Windows), this requires no compilation and merely downloads and installs a platform-specific pre-built wheel.

A Detailed Walk Through
#######################

First, save the following lines as :code:`ghz_3.py` in a folder where you want to install QMAP and run the example:

    .. code-block:: python

        from qiskit import QuantumCircuit
        from qiskit.providers.fake_provider import FakeLondon
        from mqt import qmap

        # create your quantum circuit
        circ = QuantumCircuit(3)
        circ.h(0)
        circ.cx(0, 1)
        circ.cx(0, 2)
        circ.measure_all()
        print(circ.draw(fold=-1))

        # compile circuit to 5-qubit London Architecture
        circ_mapped, results = qmap.compile(circ, FakeLondon())
        print(circ_mapped.draw(fold=-1))

        # print the result
        print("Additional Gates: %d" % (results.json()["statistics"]["additional_gates"]))

Then, the following snippet shows the installation process from setting up the virtual environment to running a small example program.

    .. code-block:: console

        $ python3 -m venv venv
        $ . venv/bin/activate
        (venv) $ pip install -U pip setuptools wheel
        (venv) $ pip install mqt.qmap
        (venv) $ python3 ghz_3.py
                ┌───┐           ░ ┌─┐
           q_0: ┤ H ├──■────■───░─┤M├──────
                └───┘┌─┴─┐  │   ░ └╥┘┌─┐
           q_1: ─────┤ X ├──┼───░──╫─┤M├───
                     └───┘┌─┴─┐ ░  ║ └╥┘┌─┐
           q_2: ──────────┤ X ├─░──╫──╫─┤M├
                          └───┘ ░  ║  ║ └╥┘
        meas: 3/═══════════════════╩══╩══╩═
                                   0  1  2
                                                            ░
                    ┌───┐┌───┐          ┌─┐
           q_0 -> 0 ┤ H ├┤ X ├──■───────┤M├──────
                    └───┘└─┬─┘┌─┴─┐     └╥┘┌─┐
           q_1 -> 1 ───────■──┤ X ├──■───╫─┤M├───
                              └───┘┌─┴─┐ ║ └╥┘┌─┐
           q_2 -> 2 ───────────────┤ X ├─╫──╫─┤M├
                                   └───┘ ║  ║ └╥┘
           q_3 -> 3 ─────────────────────╫──╫──╫─
                                         ║  ║  ║
           q_4 -> 4 ─────────────────────╫──╫──╫─
                                         ║  ║  ║
               c: 3/═════════════════════╩══╩══╩═
                                         1  0  2

        Additional Gates: 1


Building from Source for Performance
####################################

In order to get the best performance out of QMAP and enable platform-specific compiler optimizations that cannot be enabled on portable wheels, it is recommended to build the package from source via:

    .. code-block:: console

        (venv) $ pip install mqt.qmap --no-binary mqt.qmap

This requires a `C++ compiler <https://en.wikipedia.org/wiki/List_of_compilers#C++_compilers>`_ compiler supporting *C++17*, a minimum `CMake <https://cmake.org/>`_ version of *3.14* and the `SMT solver Z3 <https://github.com/Z3Prover/z3>`_. Z3 has to be installed and the dynamic linker has to be able to find the library. This can be accomplished in a multitude of ways:

- Under Ubuntu 20.04 and newer: :code:`sudo apt-get install libz3-dev`
- Under macOS: :code:`brew install z3`
- Alternatively: :code:`pip install z3-solver` and then append the corresponding path to the library path (:code:`LD_LIBRARY_PATH` under Linux, :code:`DYLD_LIBRARY_PATH` under macOS), e.g. via

    .. code-block:: console

        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(python -c "import z3; print(z3.__path__[0]+'/lib')")

- Download pre-built binaries from https://github.com/Z3Prover/z3/releases and copy the files to the respective system directories
- Build Z3 from source and install it to the system

The library is continuously tested under Linux, MacOS, and Windows using the `latest available system versions for GitHub Actions <https://github.com/actions/virtual-environments>`_.
In order to access the latest build logs, visit `qmap/actions/workflows/ci.yml <https://github.com/cda-tum/qmap/actions/workflows/ci.yml>`_.

.. note::
    We noticed some issues when compiling with Microsoft's *MSCV* compiler toolchain. If you want to start development on this project under Windows, consider using the *clang* compiler toolchain. A detailed description of how to set this up can be found `here <https://docs.microsoft.com/en-us/cpp/build/clang-support-msbuild?view=msvc-160>`_.
