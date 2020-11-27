__copyright__ = """This file is part of SCINE Utilities.
This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
See LICENSE.txt for details.
"""

from dev.conan.base import ScineConan


class ScineUtilsConan(ScineConan):
    name = "scine_utilities"
    version = "3.0.0"
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
        "microarch": ["detect", "none"]
    }
    default_options = {"shared": True, "python": False,
                       "tests": False, "coverage": False,
                       "microarch": "none"}
    exports = "dev/conan/*.py"
    exports_sources = ["dev/cmake/*", "src/*", "CMakeLists.txt", "README.rst",
                       "LICENSE.txt", "dev/conan/hook.cmake"]
    requires = ["eigen/[~=3.3.7]@conan/stable",
                "boost/[>1.65.0]@conan/stable",
                "yaml-cpp/0.6.3@scine/stable"]

    def requirements(self):
        self.requires("irc/6d5c7c37@scine/stable", private=self.options.shared)

    def _configure_cmake(self):
        return super()._configure_cmake_base("UtilsOS", None)
