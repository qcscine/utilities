/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Core/DerivedModule.h>
#include <Utils/ExternalQC/Orca/OrcaCalculator.h>
#include <Utils/Module/OSUtilsModule.h>

namespace Scine {
namespace Utils {

OSUtilsModule::OSUtilsModule() noexcept : calculator{ExternalQC::OrcaCalculator::model} {
  Core::DerivedModule::checkInvariants(*this);
}

std::string OSUtilsModule::name() const noexcept {
  return "OSUtils";
};

boost::any OSUtilsModule::get(const std::string& interface, const std::string& model) const {
  using CalculatorPtr = std::shared_ptr<Core::Calculator>;
  auto calculatorName = model;
  std::transform(calculatorName.begin(), calculatorName.end(), calculatorName.begin(), ::toupper);
  if (interface == Core::Calculator::interface) {
    if (calculatorName == ExternalQC::OrcaCalculator::model) {
      return static_cast<CalculatorPtr>(std::make_shared<ExternalQC::OrcaCalculator>());
    }
  }
  // Throw an exception if the we cannot pattern match a model
  throw Core::ClassNotImplementedError();
};

bool OSUtilsModule::has(const std::string& interface, const std::string& model) const noexcept {
  return Core::DerivedModule::has(interface, model, *this);
}

std::vector<std::string> OSUtilsModule::announceInterfaces() const noexcept {
  return Core::DerivedModule::announceInterfaces(*this);
}

std::vector<std::string> OSUtilsModule::announceModels(const std::string& interface) const noexcept {
  return Core::DerivedModule::announceModels(interface, *this);
}

std::shared_ptr<Core::Module> OSUtilsModule::make() {
  return std::make_shared<OSUtilsModule>();
}

} // namespace Utils
} // namespace Scine
