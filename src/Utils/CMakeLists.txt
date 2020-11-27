file(GLOB_RECURSE UTILS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/Utils/*.h)
file(GLOB_RECURSE UTILS_CPPS ${CMAKE_CURRENT_SOURCE_DIR}/Utils/*.cpp)

file(GLOB UTILS_MODULE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/Module/*.h)
file(GLOB UTILS_MODULE_CPPS ${CMAKE_CURRENT_SOURCE_DIR}/Module/*.cpp)

file(GLOB_RECURSE UTILS_TEST_CPPS ${CMAKE_CURRENT_SOURCE_DIR}/Tests/*.cpp)

file(GLOB_RECURSE UTILS_PYTHON_CPPS ${CMAKE_CURRENT_SOURCE_DIR}/Python/*.cpp)
file(GLOB_RECURSE UTILS_PYTHON_PYS ${CMAKE_CURRENT_SOURCE_DIR}/Python/*.py)

# Find/Import dependecies
include(ImportCore)
import_core()
include(ImportYamlCpp)
import_yamlcpp()
include(ImportIRC)
import_irc()
include(AddEigen)

if(NOT TARGET Boost::filesystem OR NOT TARGET Boost::system)
  find_package(Boost REQUIRED COMPONENTS filesystem system)
endif()

# Only look for OpenMP if parallelization is sought
if(SCINE_PARALLELIZE)
  find_package(OpenMP)
endif()

# Obey standard CMake behavior regarding shared/static libraries
add_library(UtilsOS ${UTILS_HEADERS} ${UTILS_CPPS})
if(NOT BUILD_SHARED_LIBS)
  set_target_properties(UtilsOS PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

set_target_properties(UtilsOS PROPERTIES
  OUTPUT_NAME utilsos
  WINDOWS_EXPORT_ALL_SYMBOLS ON
)

target_link_libraries(UtilsOS
  PUBLIC
    Boost::boost
    Boost::filesystem
    Boost::system
    Scine::CoreHeaders
  PRIVATE
    $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
    yaml-cpp
    irc
)

target_include_directories(UtilsOS
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include/Scine>
  PRIVATE
    ${YAML_CPP_INCLUDE_DIR}
)
add_eigen(UtilsOS PUBLIC)
add_library(Scine::UtilsOS ALIAS UtilsOS)

# Module
add_library(UtilsOSModule SHARED ${UTILS_MODULE_HEADERS} ${UTILS_MODULE_CPPS})
target_link_libraries(UtilsOSModule PUBLIC UtilsOS Scine::Core)
set_target_properties(UtilsOSModule PROPERTIES
  OUTPUT_NAME "utilsos.module"
  PREFIX ""
)

# -- Define Install
# Headers
install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Utils/
  DESTINATION include/Scine/Utils/
  FILES_MATCHING PATTERN "*.h"
)

# Targets
install(
  TARGETS UtilsOS UtilsOSModule
  EXPORT utilsOSTargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

# Component-level cmake package files
scine_install_component_cmake_files(
  COMPONENT ScineUtilsOS
  EXPORT_NAME utilsOSTargets
)

# Tests
if(SCINE_BUILD_TESTS)
  add_executable(UtilsOS_tests ${UTILS_TEST_CPPS})
  target_link_libraries(UtilsOS_tests PRIVATE
    GTest::Main
    GMock::GMock
    UtilsOS
    Scine::Core
    yaml-cpp
    ${CMAKE_DL_LIBS}
  )
  target_include_directories(UtilsOS_tests PRIVATE ${YAML_CPP_INCLUDE_DIR})
  add_test(
    NAME UtilsOSTests
    COMMAND UtilsOS_tests --gtest_output=xml:gtest_report.xml
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
  set_tests_properties(UtilsOSTests PROPERTIES
    ENVIRONMENT SCINE_MODULE_PATH=${CMAKE_CURRENT_BINARY_DIR}
  )
  # Copy files needed for tests
  file(
    COPY ${CMAKE_CURRENT_SOURCE_DIR}/Tests/Resources
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()

# Python bindings
if(SCINE_BUILD_PYTHON_BINDINGS)
  include(ImportPybind11)
  import_pybind11()
  set(PYBIND11_PYTHON_VERSION ${PYTHONVERSION})

  # Python module
  pybind11_add_module(scine_utilities ${UTILS_PYTHON_CPPS})
  set_target_properties(scine_utilities PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/scine_utilities
  )
  target_link_libraries(scine_utilities PRIVATE Scine::UtilsOS Scine::Core yaml-cpp)
  target_include_directories(scine_utilities PRIVATE ${YAML_CPP_INCLUDE_DIR})
  install(CODE
    "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install --prefix=${CMAKE_INSTALL_PREFIX} --upgrade ${CMAKE_CURRENT_BINARY_DIR}
                     RESULT_VARIABLE retcode)
     if(NOT \$\{retcode\} EQUAL 0)
       message(FATAL_ERROR \"Fatal error when installing Python module using PIP.\")
     endif()"
  )

  # Python tests
  if(SCINE_BUILD_TESTS)
    add_test(
      NAME UtilsOS_Python_Bindings
      COMMAND ${PYTHON_EXECUTABLE} -B -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/Python/Tests --junitxml=${CMAKE_CURRENT_BINARY_DIR}/pytest_report.xml
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(UtilsOS_Python_Bindings PROPERTIES
      ENVIRONMENT PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}
    )
  endif()

  # Copy the setup file and README.rst into the binary dir so you can install
  # using pip or make a distribution with the compiled binary
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/Python/pkginit.py
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/scine_utilities
  )
  file(RENAME
    ${CMAKE_CURRENT_BINARY_DIR}/scine_utilities/pkginit.py
    ${CMAKE_CURRENT_BINARY_DIR}/scine_utilities/__init__.py
  )
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/Python/README.rst DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  set(LIBNAME_CORE "'${CMAKE_SHARED_LIBRARY_PREFIX}core${CMAKE_SHARED_LIBRARY_SUFFIX}'")
  add_custom_command(TARGET UtilsOSModule POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:UtilsOSModule> ${CMAKE_CURRENT_BINARY_DIR}/scine_utilities
    COMMENT "Copying UtilsOS module into python package directory"
  )
  if(BUILD_SHARED_LIBS)
    add_custom_command(TARGET UtilsOS POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:UtilsOS> $<TARGET_FILE:Scine::Core> ${CMAKE_CURRENT_BINARY_DIR}/scine_utilities
      COMMENT "Copying dependent shared libraries into python package directory"
    )
    set(_utils_name "${CMAKE_SHARED_LIBRARY_PREFIX}utilsos${CMAKE_SHARED_LIBRARY_SUFFIX}")
    include(TargetLibName)
    target_lib_filename(Scine::Core _core_name)
    message(STATUS "deduced library names: ${_utils_name}, ${_core_name}")
    set(utils_PY_DEPS ", \"${_utils_name}\", \"${_core_name}\"")
    unset(_utils_name)
    unset(_core_name)
  endif()

  # Add setuptools file
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Python/setup.py
    ${CMAKE_CURRENT_BINARY_DIR}/setup.py
    @ONLY
  )
  unset(utils_PY_DEPS)
  file(
    COPY ${PROJECT_SOURCE_DIR}/README.rst
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  )

  # Python documentation
  include(SphinxDocumentation)
  scine_sphinx_documentation(
    TARGET scine_utilities
    CONFIGURATION ${CMAKE_CURRENT_SOURCE_DIR}/Python/Sphinx/conf.py.in
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Python/Sphinx
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/doc-py
    DOCTEST
  )
  # The UtilsOSModule is technically also a dependency of the documentation
  if(TARGET scine_utilitiesDocumentation)
    add_dependencies(scine_utilitiesDocumentation UtilsOSModule)
  endif()
endif()

if(SCINE_BUILD_PYTHON_BINDINGS)
  file(COPY ${PROJECT_SOURCE_DIR}/LICENSE.txt  ${EXTERNAL_LICENSE_FILES}
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/scine_utilities
  )
  file(COPY ${EXTERNAL_LICENSE_FILES}
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/external_licenses
  )
endif()

if(WIN32)
  set_target_properties(UtilsOS PROPERTIES PREFIX "" IMPORT_PREFIX "")
  if(MINGW)
    target_link_libraries(UtilsOS PRIVATE ws2_32)
  endif()
endif()
