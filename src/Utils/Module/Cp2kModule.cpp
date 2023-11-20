/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Cp2kModule.h"
#include <Core/DerivedModule.h>
#include <Utils/ExternalQC/Cp2k/Cp2kCalculator.h>

namespace Scine {
namespace Utils {

using Cp2kInterfaceModelMap =
    boost::mpl::map<boost::mpl::pair<Scine::Core::Calculator, boost::mpl::vector<ExternalQC::Cp2kCalculator>>>;

std::string Cp2kModule::name() const noexcept {
  return "Cp2k";
}

bool Cp2kModule::cp2kFound() {
  return std::getenv("CP2K_BINARY_PATH") != nullptr;
}

boost::any Cp2kModule::get(const std::string& interface, const std::string& model) const {
  if (interface == Scine::Core::Calculator::interface && model == ExternalQC::Cp2kCalculator::model && !cp2kFound()) {
    throw Scine::Core::ClassNotImplementedError();
  }

  boost::any resolved = Scine::Core::DerivedModule::resolve<Cp2kInterfaceModelMap>(interface, model);
  // Throw an exception if we could not match an interface or model
  if (resolved.empty()) {
    throw Scine::Core::ClassNotImplementedError();
  }
  return resolved;
}

bool Cp2kModule::has(const std::string& interface, const std::string& model) const noexcept {
  if (interface == Scine::Core::Calculator::interface) {
    if (model == ExternalQC::Cp2kCalculator::model) {
      return cp2kFound();
    }
  }

  return Scine::Core::DerivedModule::has<Cp2kInterfaceModelMap>(interface, model);
}

std::vector<std::string> Cp2kModule::announceInterfaces() const noexcept {
  return Scine::Core::DerivedModule::announceInterfaces<Cp2kInterfaceModelMap>();
}

std::vector<std::string> Cp2kModule::announceModels(const std::string& interface) const noexcept {
  auto models = Scine::Core::DerivedModule::announceModels<Cp2kInterfaceModelMap>(interface);

  auto cp2kModel = ExternalQC::Cp2kCalculator::model;
  if (interface == Scine::Core::Calculator::interface && !cp2kFound()) {
    models.erase(std::remove(std::begin(models), std::end(models), cp2kModel), std::end(models));
  }

  return models;
}

std::shared_ptr<Scine::Core::Module> Cp2kModule::make() {
  return std::make_shared<Cp2kModule>();
}

} // namespace Utils
} // namespace Scine
