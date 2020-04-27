file(GLOB_RECURSE UTILS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/Utils/*.h)
file(GLOB_RECURSE UTILS_CPPS ${CMAKE_CURRENT_SOURCE_DIR}/Utils/*.cpp)
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
set(Boost_USE_STATIC_LIBS OFF)
find_package(Boost REQUIRED COMPONENTS filesystem log system)
find_package(OpenMP)

# Shared library
add_library(UtilsOS SHARED ${UTILS_CPPS} ${UTILS_HEADERS})
set_target_properties(UtilsOS PROPERTIES
  OUTPUT_NAME utilsos
  WINDOWS_EXPORT_ALL_SYMBOLS ON
  POSITION_INDEPENDENT_CODE ON
)
target_compile_options(UtilsOS PRIVATE $<$<BOOL:${OpenMP_CXX_FOUND}>:${OpenMP_CXX_FLAGS}>)
# Include Core only as there should be no link dependency on
#  Scine::Core for this target.
target_include_directories(UtilsOS
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
    $<TARGET_PROPERTY:Scine::Core,INTERFACE_INCLUDE_DIRECTORIES>
    $<TARGET_PROPERTY:yaml-cpp,INCLUDE_DIRECTORIES>
  PRIVATE
    ${LIBIRC_INCLUDE_DIR}
)
target_compile_definitions(UtilsOS PUBLIC BOOST_ALL_DYN_LINK)
target_link_libraries(UtilsOS
  PUBLIC
    Boost::boost
    Boost::filesystem
    Boost::log
    Boost::system
    Scine::Core
    yaml-cpp
  PRIVATE
    $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
)
add_eigen(UtilsOS PUBLIC)
target_compile_options(UtilsOS
  PUBLIC
    $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>
)



# Add namespaced aliases
add_library(Scine::UtilsOS ALIAS UtilsOS)

# -- Define Install
# Headers
install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Utils/
  DESTINATION include/Utils
  FILES_MATCHING PATTERN "*.h"
)

# Targets
install(
  TARGETS UtilsOS
  EXPORT utilsOSTargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include/Utils
)

# Symlink on install
install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_INSTALL_PREFIX}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}utilsos${CMAKE_SHARED_LIBRARY_SUFFIX} \
        ${CMAKE_INSTALL_PREFIX}/lib/utilsos.module${CMAKE_SHARED_LIBRARY_SUFFIX})"
)

# Symlink on build
add_custom_command(TARGET UtilsOS POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E create_symlink
                   "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}utilsos${CMAKE_SHARED_LIBRARY_SUFFIX}"
                   "${CMAKE_CURRENT_BINARY_DIR}/utilsos.module${CMAKE_SHARED_LIBRARY_SUFFIX}"
)

# Component-level cmake package files
scine_install_component_cmake_files(
  COMPONENT UtilsOS
  EXPORT_NAME utilsOSTargets
)

# Tests
if(SCINE_BUILD_TESTS)
  add_executable(UtilsOS_tests ${UTILS_TEST_CPPS})
  target_link_libraries(UtilsOS_tests
    PUBLIC
      gtest_main
      gmock
      ${PROJECT_NAME}
  )
  target_compile_options(UtilsOS_tests
    PUBLIC
      $<TARGET_PROPERTY:UtilsOS,INTERFACE_COMPILE_OPTIONS>
  )
  add_test(
    NAME UtilsOSTests
    COMMAND UtilsOS_tests --gtest_output=xml:gtest_report.xml
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
  set_tests_properties(UtilsOSTests PROPERTIES
    ENVIRONMENT SCINE_MODULE_PATH=${CMAKE_CURRENT_BINARY_DIR}
  )
  # Hard-code location of the ORCA output files into the test by configuring a header file for it.
  set(externalQCTests_dir "${CMAKE_CURRENT_SOURCE_DIR}/Tests/ExternalQC/")
  configure_file(${externalQCTests_dir}/externalQC_output_location.h.in
                 ${externalQCTests_dir}/externalQC_output_location.h)

  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/Tests/IO/traj.xyz
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/test_resources/)
endif()

# Python bindings
if(SCINE_BUILD_PYTHON_BINDINGS)
  include(ImportPybind11)
  import_pybind11()
  set(PYBIND11_PYTHON_VERSION ${PYTHONVERSION})

  # Python module
  pybind11_add_module(scine_utilities ${UTILS_PYTHON_CPPS})
  set_target_properties(scine_utilities PROPERTIES SUFFIX ".so")
  set_target_properties(scine_utilities PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/scine_utilities)
  target_link_libraries(scine_utilities PRIVATE Scine::UtilsOS Scine::Core)
  target_compile_options(scine_utilities
    PUBLIC
      $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>
  )
  install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install --prefix=${CMAKE_INSTALL_PREFIX} --upgrade ${CMAKE_CURRENT_BINARY_DIR})")

  # Python tests
  if(SCINE_BUILD_TESTS)
    add_test(
      NAME UtilsOS_Python_Bindings
      COMMAND ${PYTHON_EXECUTABLE} -B -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/Python/Tests --junitxml=${CMAKE_CURRENT_BINARY_DIR}/pytest_report.xml
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(UtilsOS_Python_Bindings PROPERTIES
      ENVIRONMENT SCINE_MODULE_PATH=${CMAKE_CURRENT_BINARY_DIR}
      ENVIRONMENT PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}
    )
  endif()

  # Copy the setup file and README.md into the binary dir so you can install
  # using pip or make a distribution with the compiled binary
  file(
    COPY
      ${CMAKE_CURRENT_SOURCE_DIR}/Python/__init__.py
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/scine_utilities
  )
  file(
    COPY
      ${CMAKE_CURRENT_SOURCE_DIR}/Python/setup.py
      ${PROJECT_SOURCE_DIR}/README.md
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()

if(WIN32)
  set_target_properties(${PROJECT_NAME} PROPERTIES PREFIX "")
  set_target_properties(${PROJECT_NAME} PROPERTIES IMPORT_PREFIX "")
  set_target_properties(${PROJECT_NAME} PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

# Add link library specific for MinGW
if(WIN32)
  target_link_libraries(UtilsOS PRIVATE $<$<BOOL:MINGW>:ws2_32>)
endif()
