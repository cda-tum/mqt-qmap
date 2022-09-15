import os
import re
import subprocess
import sys
from contextlib import suppress
from pathlib import Path

from setuptools import Extension, find_namespace_packages, setup
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        super().__init__(name, sources=[])
        self.sourcedir = str(Path(sourcedir).resolve())


class CMakeBuild(build_ext):
    def build_extension(self, ext):
        from setuptools_scm import get_version

        version = get_version(root=".", relative_to=__file__)

        extdir = str(Path(self.get_ext_fullpath(ext.name)).parent.resolve())

        cmake_generator = os.environ.get("CMAKE_GENERATOR", "")
        cfg = "Debug" if self.debug else "Release"

        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DQMAP_VERSION_INFO={version}",
            f"-DCMAKE_BUILD_TYPE={cfg}",
            "-DBINDINGS=ON",
        ]
        build_args = []

        if self.compiler.compiler_type != "msvc":
            if not cmake_generator:
                cmake_args += ["-GNinja"]
        else:
            # Single config generators are handled "normally"
            single_config = any(x in cmake_generator for x in {"NMake", "Ninja"})
            # CMake allows an arch-in-generator style for backward compatibility
            contains_arch = any(x in cmake_generator for x in {"ARM", "Win64"})
            # Convert distutils Windows platform specifiers to CMake -A arguments
            plat_to_cmake = {
                "win32": "Win32",
                "win-amd64": "x64",
                "win-arm32": "ARM",
                "win-arm64": "ARM64",
            }
            # Specify the arch if using MSVC generator, but only if it doesn't
            # contain a backward-compatibility arch spec already in the
            # generator name.
            if not single_config and not contains_arch:
                cmake_args += ["-A", plat_to_cmake[self.plat_name]]
            # Multi-config generators have a different way to specify configs
            if not single_config:
                cmake_args += [f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}"]
                build_args += ["--config", cfg]

        # cross-compile support for macOS - respect ARCHFLAGS if set
        if sys.platform.startswith("darwin"):
            archs = re.findall(r"-arch (\S+)", os.environ.get("ARCHFLAGS", ""))
            if archs:
                cmake_args += ["-DCMAKE_OSX_ARCHITECTURES={}".format(";".join(archs))]

        # Set CMAKE_BUILD_PARALLEL_LEVEL to control the parallel build level across all generators.
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ and hasattr(self, "parallel") and self.parallel:
            build_args += [f"-j{self.parallel}"]

        if sys.platform == "win32":
            cmake_args += ["-T", "ClangCl"]

        build_dir = Path(self.build_temp)
        build_dir.mkdir(parents=True, exist_ok=True)
        with suppress(FileNotFoundError):
            Path(build_dir / "CMakeCache.txt").unlink()

        subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp)
        subprocess.check_call(
            ["cmake", "--build", ".", "--target", ext.name.split(".")[-1]] + build_args,
            cwd=self.build_temp,
        )


README_PATH = Path(__file__).parent.resolve() / "README.md"
with README_PATH.open() as readme_file:
    README = readme_file.read()

setup(
    name="mqt.qmap",
    author="Lukas Burgholzer",
    author_email="lukas.burgholzer@jku.at",
    description="A tool for Quantum Circuit Mapping",
    long_description=README,
    long_description_content_type="text/markdown",
    python_requires=">=3.7",
    license="MIT",
    url="https://www.cda.cit.tum.de/research/ibm_qx_mapping/",
    ext_modules=[CMakeExtension("mqt.qmap.pyqmap")],
    cmdclass={"build_ext": CMakeBuild},
    zip_safe=False,
    packages=find_namespace_packages(include=["mqt.*"]),
    install_requires=["qiskit-terra>=0.20.2,<0.22.0", "retworkx>=0.11.0", "networkx>=1.9"],
    extras_require={
        "test": ["pytest~=7.1.1", "mqt.qcec~=2.0.0rc7"],
        "coverage": ["coverage[toml]~=6.4.2", "pytest-cov~=3.0.0"],
        "docs": [
            "sphinx>=5.1.1",
            "sphinx-rtd-theme",
            "sphinxcontrib-bibtex~=2.5",
            "sphinx-copybutton",
            "sphinx-hoverxref~=1.1.3",
            "pybtex>=0.24",
            "importlib_metadata>=3.6; python_version < '3.10'",
        ],
        "dev": ["mqt.qmap[test, coverage, docs]"],  # requires Pip 21.2 or newer
    },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: C++",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: MacOS",
        "Operating System :: POSIX :: Linux",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Electronic Design Automation (EDA)",
    ],
    keywords="MQT quantum compilation mapping",
    project_urls={
        "Source": "https://github.com/cda-tum/qmap/",
        "Tracker": "https://github.com/cda-tum/qmap/issues",
        "Research": "https://www.cda.cit.tum.de/research/ibm_qx_mapping/",
        "Documentation": "https://mqtqmap.readthedocs.io",
    },
)
