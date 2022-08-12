__copyright__ = """This file is part of SCINE Utilities.
This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
See LICENSE.txt for details.
"""

import sys
from dev.conan.base import ScineConan

class ScineUtilsConan(ScineConan):
    name = "scine_utilities"
    version = "6.0.0"
    url = "https://github.com/qcscine/utilities"
    description = """
Functionality which is used in most SCINE modules. It is vital for the correct
functioning of all SCINE modules, but it does not directly provide any
functionality for end users. Therefore, only developers of SCINE need to
directly interact with this repository."""
    options = {
        "shared": [True, False],
        "python": [True, False],
        "tests": [True, False],
        "coverage": [True, False],
        "microarch": ["detect", "none"],
        "python_version": "ANY"
    }
    python_version_string = str(sys.version_info.major) + \
        "." + str(sys.version_info.minor) 
    default_options = {"shared": True, "python": False,
                       "tests": False, "coverage": False,
                       "microarch": "none", "python_version": python_version_string}
    exports = "dev/conan/*.py"
    exports_sources = ["dev/cmake/*", "src/*", "CMakeLists.txt", "README.rst",
                       "LICENSE.txt", "dev/conan/hook.cmake",
                       "dev/conan/glue/*"]
    requires = ["eigen/[~=3.3.7]", "boost/[>1.65.0]", "yaml-cpp/0.6.3", "lbfgspp/0.1.0"]
    cmake_name = "UtilsOS"
    cmake_definitions = {
        "CMAKE_UNITY_BUILD": "ON",
        "CMAKE_UNITY_BUILD_BATCH_SIZE": 16
    }

    def requirements(self):
        # Cannot mark irc as private until #7585 gets merged and a new conan
        # was released. Replace with:
        # >>>  self.requires("irc/6d5c7c37", private=self.options.shared)
        self.requires("irc/6d5c7c37")

    def package_info(self):
        super().package_info()

        # Add Core's targets as components and model dependency chain
        self.cpp_info.components["CoreHeaders"].includedirs = ["include/Scine"]
        self.cpp_info.components["CoreHeaders"].requires = ["boost::boost"]
        self.cpp_info.components["Core"].libs = ["core"]
        self.cpp_info.components["Core"].requires = ["CoreHeaders"]
        self.cpp_info.components["UtilsOS"].libs.remove("core")
        self.cpp_info.components["UtilsOS"].requires.append("CoreHeaders")
        self.cpp_info.components["UtilsOS"].cxxflags = ["-fopenmp"]
        self.cpp_info.components["UtilsOS"].sharedlinkflags = ["-fopenmp"]
        self.cpp_info.components["UtilsOS"].exelinkflags = ["-fopenmp"]
