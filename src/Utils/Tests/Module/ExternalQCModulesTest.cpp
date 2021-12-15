/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Core/Interfaces/Calculator.h>
#include <Core/ModuleManager.h>
#include <Utils/ExternalQC/Cp2k/Cp2kCalculator.h>
#include <Utils/ExternalQC/ExternalProgram.h>
#include <Utils/ExternalQC/Gaussian/GaussianCalculator.h>
#include <Utils/ExternalQC/Orca/OrcaCalculator.h>
#include <Utils/ExternalQC/Turbomole/TurbomoleCalculator.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

class ExternalQCModulesTest : public Test {};

TEST_F(ExternalQCModulesTest, GaussianModuleIsWorkingCorrectly) {
#ifndef _WIN32
  const char* envVariablePtr = std::getenv("GAUSSIAN_BINARY_PATH");
  if (envVariablePtr) {
    unsetenv("GAUSSIAN_BINARY_PATH"); // only unset the variable if it was previously set
  }

  auto& manager = Core::ModuleManager::getInstance();
  auto modelName = ExternalQC::GaussianCalculator::model;

  auto loadedModels = manager.getLoadedModels(Core::Calculator::interface);
  bool gaussianAvailable = std::find(loadedModels.begin(), loadedModels.end(), modelName) != loadedModels.end();
  ASSERT_FALSE(gaussianAvailable);
  ASSERT_FALSE(manager.has(Core::Calculator::interface, modelName));
  EXPECT_THROW(manager.get<Core::Calculator>(modelName, "Gaussian"), Core::ClassNotImplementedError);
  EXPECT_THROW(manager.get<Core::Calculator>(Core::Calculator::supports("DFT"), "Gaussian"), std::runtime_error);

  setenv("GAUSSIAN_BINARY_PATH", "g09", true);

  loadedModels = manager.getLoadedModels(Core::Calculator::interface);
  gaussianAvailable = std::find(loadedModels.begin(), loadedModels.end(), modelName) != loadedModels.end();
  ASSERT_TRUE(gaussianAvailable);
  ASSERT_TRUE(manager.has(Core::Calculator::interface, modelName));

  // These two lines should not throw an exception
  auto calc1 = manager.get<Core::Calculator>(modelName, "Gaussian");
  auto calc2 = manager.get<Core::Calculator>(Core::Calculator::supports("DFT"), "Gaussian");

  // Clean up
  if (!envVariablePtr) {
    unsetenv("GAUSSIAN_BINARY_PATH"); // only unset the variable if it was set to g09
  }
  else {
    setenv("GAUSSIAN_BINARY_PATH", envVariablePtr, true);
  }
#endif
}

TEST_F(ExternalQCModulesTest, OrcaModuleIsWorkingCorrectly) {
#ifndef _WIN32
  const char* envVariablePtr = std::getenv("ORCA_BINARY_PATH");
  if (envVariablePtr) {
    unsetenv("ORCA_BINARY_PATH");
  }

  auto& manager = Core::ModuleManager::getInstance();
  auto modelName = ExternalQC::OrcaCalculator::model;

  auto loadedModels = manager.getLoadedModels(Core::Calculator::interface);
  bool orcaAvailable = std::find(loadedModels.begin(), loadedModels.end(), modelName) != loadedModels.end();
  ASSERT_FALSE(orcaAvailable);
  ASSERT_FALSE(manager.has(Core::Calculator::interface, modelName));
  EXPECT_THROW(manager.get<Core::Calculator>(modelName, "Orca"), Core::ClassNotImplementedError);
  EXPECT_THROW(manager.get<Core::Calculator>(Core::Calculator::supports("DFT"), "Orca"), std::runtime_error);

  setenv("ORCA_BINARY_PATH", "orca", true);

  loadedModels = manager.getLoadedModels(Core::Calculator::interface);
  orcaAvailable = std::find(loadedModels.begin(), loadedModels.end(), modelName) != loadedModels.end();
  ASSERT_TRUE(orcaAvailable);
  ASSERT_TRUE(manager.has(Core::Calculator::interface, modelName));

  // These two lines should not throw an exception
  auto calc1 = manager.get<Core::Calculator>(modelName, "Orca");
  auto calc2 = manager.get<Core::Calculator>(Core::Calculator::supports("DFT"), "Orca");

  // Clean up
  if (!envVariablePtr) {
    unsetenv("ORCA_BINARY_PATH"); // only unset the variable if it was set to g09
  }
  else {
    setenv("ORCA_BINARY_PATH", envVariablePtr, true);
  }
#endif
}

TEST_F(ExternalQCModulesTest, Cp2kModuleIsWorkingCorrectly) {
#ifndef _WIN32
  const char* envVariablePtr = std::getenv("CP2K_BINARY_PATH");
  if (envVariablePtr) {
    unsetenv("CP2K_BINARY_PATH");
  }

  auto& manager = Core::ModuleManager::getInstance();
  auto modelName = ExternalQC::Cp2kCalculator::model;

  auto loadedModels = manager.getLoadedModels(Core::Calculator::interface);
  bool cp2kAvailable = std::find(loadedModels.begin(), loadedModels.end(), modelName) != loadedModels.end();
  ASSERT_FALSE(cp2kAvailable);
  ASSERT_FALSE(manager.has(Core::Calculator::interface, modelName));
  EXPECT_THROW(manager.get<Core::Calculator>(modelName, "Cp2k"), Core::ClassNotImplementedError);
  EXPECT_THROW(manager.get<Core::Calculator>(Core::Calculator::supports("DFT"), "Cp2k"), std::runtime_error);

  setenv("CP2K_BINARY_PATH", "cp2k", true);

  loadedModels = manager.getLoadedModels(Core::Calculator::interface);
  cp2kAvailable = std::find(loadedModels.begin(), loadedModels.end(), modelName) != loadedModels.end();
  ASSERT_TRUE(cp2kAvailable);
  ASSERT_TRUE(manager.has(Core::Calculator::interface, modelName));

  // These two lines should not throw an exception
  auto calc1 = manager.get<Core::Calculator>(modelName, "Cp2k");
  auto calc2 = manager.get<Core::Calculator>(Core::Calculator::supports("DFT"), "Cp2k");

  // Clean up
  if (!envVariablePtr) {
    unsetenv("CP2K_BINARY_PATH");
  }
  else {
    setenv("CP2K_BINARY_PATH", envVariablePtr, true);
  }
#endif
}

TEST_F(ExternalQCModulesTest, TurbomoleModuleIsWorkingCorrectly) {
#ifndef _WIN32
  const char* envVariablePtr = std::getenv("TURBODIR");
  if (envVariablePtr) {
    unsetenv("TURBODIR");
  }

  auto& manager = Core::ModuleManager::getInstance();
  auto modelName = ExternalQC::TurbomoleCalculator::model;

  auto loadedModels = manager.getLoadedModels(Core::Calculator::interface);
  bool turbomoleAvailable = std::find(loadedModels.begin(), loadedModels.end(), modelName) != loadedModels.end();
  ASSERT_FALSE(turbomoleAvailable);
  ASSERT_FALSE(manager.has(Core::Calculator::interface, modelName));

  EXPECT_THROW(manager.get<Core::Calculator>(modelName, "Turbomole"), Core::ClassNotImplementedError);
  EXPECT_THROW(manager.get<Core::Calculator>(Core::Calculator::supports("DFT"), "Turbomole"), std::runtime_error);

  setenv("TURBODIR", "dscf", true);

  loadedModels = manager.getLoadedModels(Core::Calculator::interface);
  turbomoleAvailable = std::find(loadedModels.begin(), loadedModels.end(), modelName) != loadedModels.end();
  ASSERT_TRUE(turbomoleAvailable);
  ASSERT_TRUE(manager.has(Core::Calculator::interface, modelName));

  // These two lines should not throw an exception
  auto calc1 = manager.get<Core::Calculator>(modelName, "Turbomole");
  auto calc2 = manager.get<Core::Calculator>(Core::Calculator::supports("DFT"), "Turbomole");

  // Clean up
  if (!envVariablePtr) {
    unsetenv("TURBODIR");
  }
  else {
    setenv("TURBODIR", envVariablePtr, true);
  }
#endif
}

TEST_F(ExternalQCModulesTest, WorkingDirectoryCanBeSet) {
  std::string workingDir = "/home/path/to/testDir";
  ExternalQC::ExternalProgram externalProgram;
  externalProgram.setWorkingDirectory(workingDir);
  ASSERT_THAT(externalProgram.getWorkingDirectory(), NativeFilenames::addTrailingSeparator(workingDir));
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
