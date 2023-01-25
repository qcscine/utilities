__copyright__ = """This file is part of SCINE Utilities.
This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
See LICENSE.txt for details.
"""

from conans import ConanFile, CMake
import sys


class TestPackageConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "cmake_find_package"
    # CMake 3.20.0 has a bug when finding/linking openmp
    #  so it is exluded from the allowed versions, the version
    #  gap can be increased if 3.20.1 does not fix the error
    build_requires = "cmake/[>=3.18.0 <3.20.0 || >3.20.0]"
    exports_sources = "CMakeLists.txt", "test.cpp"

    def _configure(self):
        cmake = CMake(self, cmake_program='/usr/bin/cmake')
        cmake.configure()
        return cmake

    def build(self):
        cmake = self._configure()
        cmake.build()

    def test(self):
        cmake = self._configure()
        try:
            cmake.test(output_on_failure=True)
        except Exception as e:
            sys.stderr.write(str(e))
            raise e

        if self.options["scine_utilities"].python:
            self.output.info("Trying to import 'scine_utilities'")
            try:
                import scine_utilities
            except Exception as e:
                sys.stderr.write(str(e))
                raise e
            self.output.info("Import worked!")
