cmake_minimum_required(VERSION 3.9)

# Set the project name. The chosen project name must be in CamelCase. Your CMake
# tree must then provide a properly namespaced target with the same name as
# your project.
project(UtilsOS
  VERSION 5.0.0
  DESCRIPTION "Utilities to be used (statically linked) in all other SCINE modules."
)

# Set the module path for universal cmake files inclusion
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/dev/cmake)

# Component setup
include(ComponentSetup)
scine_setup_component()

# Enable testing
if(SCINE_BUILD_TESTS)
  enable_testing()
endif()

# Set rpaths of shared libraries if generated
if(APPLE)
  set(CMAKE_BUILD_RPATH "@loader_path")
  set(CMAKE_INSTALL_RPATH "@loader_path/../lib;@loader_path")
elseif(UNIX)
  set(CMAKE_BUILD_RPATH "\$ORIGIN")
  set(CMAKE_INSTALL_RPATH "\$ORIGIN/../lib;\$ORIGIN")
endif()

# Subdirectories
add_subdirectory(src)

# Documentation
scine_component_documentation(CoreDocumentation)
