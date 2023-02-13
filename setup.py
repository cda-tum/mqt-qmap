"""Setup script for the MQT QCEC package."""

import os
import re
import subprocess
import sys
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    """Class that wraps a CMake extension."""

    def __init__(self, name: str, sourcedir: str = "") -> None:
        """Initialize the CMake extension.

        Args:
        ----
        name: The name of the extension.
        sourcedir: The path to the source directory.
        """
        super().__init__(name, sources=[])
        self.sourcedir = str(Path(sourcedir).resolve())


class CMakeBuild(build_ext):
    """Class that builds a CMake extension."""

    def build_extension(self, ext: CMakeExtension) -> None:
        """Build the CMake extension.

        Args:
        ----
        ext: The CMake extension to build.
        """
        from setuptools_scm import get_version  # type: ignore[import]

        version = get_version(root=".", relative_to=__file__)

        extdir = str(Path(self.get_ext_fullpath(ext.name)).parent.resolve())  # type: ignore[no-untyped-call]

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
        Path(build_dir / "CMakeCache.txt").unlink(missing_ok=True)

        subprocess.check_call(["cmake", ext.sourcedir, *cmake_args], cwd=self.build_temp)
        subprocess.check_call(
            ["cmake", "--build", ".", "--target", ext.name.split(".")[-1], *build_args],
            cwd=self.build_temp,
        )


setup(
    ext_modules=[CMakeExtension("mqt.qmap.pyqmap")],
    cmdclass={"build_ext": CMakeBuild},
)
