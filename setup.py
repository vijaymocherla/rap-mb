#!/usr/bin/env python
"""rap_mb setup
"""
from setuptools import find_packages

import os
import sys
import platform
import subprocess

from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    build_ext.user_options = build_ext.user_options + [
        # Notes: the first option is the option string
        #        the second option is an abbreviated form of an option, which we avoid with None
        ("code-coverage", None, "enable code coverage")
    ]

    def initialize_options(self):
        self.code_coverage = "OFF"
        return build_ext.initialize_options(self)

    def run(self):
        try:
            out = subprocess.check_output(["cmake", "--version"])
            print(f"CMake version: {out}")
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: "
                + ", ".join(e.name for e in self.extensions)
            )

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.join(
            os.path.dirname(self.get_ext_fullpath(ext.name)), "rap_mb"
        )
        extdir = os.path.abspath(extdir)

        print(f"    CODE_COVERAGE = {str(self.code_coverage).upper()}")

        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + extdir,
            "-DPYTHON_EXECUTABLE=" + sys.executable,
        ]

        cfg = "Debug" if self.debug else "Release"
        build_args = ["--config", cfg]

        if platform.system() == "Windows":
            cmake_args += [
                "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}".format(cfg.upper(), extdir)
            ]
            if sys.maxsize > 2 ** 32:
                cmake_args += ["-A", "x64"]
            build_args += ["--", "/m"]
        else:
            cmake_args += ["-DCMAKE_BUILD_TYPE=" + cfg]
            build_args += ["--", "-j2"]

        env = os.environ.copy()
        env["CXXFLAGS"] = "{} -DVERSION_INFO=\\'{}\\'".format(
            env.get("CXXFLAGS", ""), self.distribution.get_version())
        cmake_args += [f"-DCODE_COVERAGE={str(self.code_coverage).upper()}"]
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        subprocess.check_call(["cmake"] + cmake_args)
        subprocess.check_call(["cmake", "--build", ".", "-j2"] + build_args)

        print()  # Add empty line for nicer output


setup(
    name="rap_mb",
    version="0.1.0",
    author="Sai Vijay Mocherla",
    description="",
    long_description="",
    packages=find_packages(exclude=['tests']),
    package_dir={"": "."},
    # add an extension module named 'rap_mb' to the package 'rap_mb'
    ext_modules=[CMakeExtension("rap_mb")],
    # add custom build_ext command
    cmdclass=dict(build_ext=CMakeBuild),
    test_suite="tests",
    zip_safe=False,
)

