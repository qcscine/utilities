/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "TestModule.h"
#include <Core/DerivedModule.h>
#include <Utils/CalculatorBasics/TestCalculator.h>

namespace Scine {
namespace Utils {

using InterfaceModelMap = boost::mpl::map<boost::mpl::pair<Scine::Core::Calculator, boost::mpl::vector<TestCalculator>>>;

std::string TestModule::name() const noexcept {
  return "Test";
}

boost::any TestModule::get(const std::string& interface, const std::string& model) const {
  boost::any resolved = Scine::Core::DerivedModule::resolve<InterfaceModelMap>(interface, model);
  // Throw an exception if we could not match an interface or model
  if (resolved.empty()) {
    throw Scine::Core::ClassNotImplementedError();
  }
  return resolved;
}

bool TestModule::has(const std::string& interface, const std::string& model) const noexcept {
  if (interface == Scine::Core::Calculator::interface) {
    if (model == TestCalculator::model) {
      return true;
    }
  }

  return Scine::Core::DerivedModule::has<InterfaceModelMap>(interface, model);
}

std::vector<std::string> TestModule::announceInterfaces() const noexcept {
  return Scine::Core::DerivedModule::announceInterfaces<InterfaceModelMap>();
}

std::vector<std::string> TestModule::announceModels(const std::string& interface) const noexcept {
  return Scine::Core::DerivedModule::announceModels<InterfaceModelMap>(interface);
}

std::shared_ptr<Scine::Core::Module> TestModule::make() {
  return std::make_shared<TestModule>();
}

} // namespace Utils
} // namespace Scine
