__copyright__ = """This file is part of SCINE Utilities.
This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
See LICENSE.txt for details.
"""

import os
from conans import ConanFile, CMake


def conan_paths_str(build_folder):
    return os.path.join(build_folder, "conan_paths.cmake")


class TestPackageConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "cmake_paths"
    build_requires = [("cmake/[>3.13.4]@scine/stable")]
    exports_sources = "CMakeLists.txt", "test.cpp"

    def _configure(self):
        cmake = CMake(self)
        cmake.definitions["CMAKE_PROJECT_UtilsTestPackage_INCLUDE"] = conan_paths_str(
            self.build_folder)
        cmake.configure()
        return cmake

    def build(self):
        cmake = self._configure()
        cmake.build()

    def test(self):
        cmake = self._configure()
        cmake.test()
