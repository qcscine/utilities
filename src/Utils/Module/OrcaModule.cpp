/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "OrcaModule.h"
#include <Core/DerivedModule.h>
#include <Utils/ExternalQC/Orca/OrcaCalculator.h>

namespace Scine {
namespace Utils {

using OrcaInterfaceModelMap =
    boost::mpl::map<boost::mpl::pair<Scine::Core::Calculator, boost::mpl::vector<ExternalQC::OrcaCalculator>>>;

std::string OrcaModule::name() const noexcept {
  return "Orca";
}

bool OrcaModule::orcaFound() {
  return std::getenv("ORCA_BINARY_PATH") != nullptr;
}

boost::any OrcaModule::get(const std::string& interface, const std::string& model) const {
  if (interface == Scine::Core::Calculator::interface && model == ExternalQC::OrcaCalculator::model && !orcaFound()) {
    throw Scine::Core::ClassNotImplementedError();
  }

  boost::any resolved = Scine::Core::DerivedModule::resolve<OrcaInterfaceModelMap>(interface, model);
  // Throw an exception if we could not match an interface or model
  if (resolved.empty()) {
    throw Scine::Core::ClassNotImplementedError();
  }
  return resolved;
}

bool OrcaModule::has(const std::string& interface, const std::string& model) const noexcept {
  if (interface == Scine::Core::Calculator::interface) {
    if (model == ExternalQC::OrcaCalculator::model) {
      return orcaFound();
    }
  }

  return Scine::Core::DerivedModule::has<OrcaInterfaceModelMap>(interface, model);
}

std::vector<std::string> OrcaModule::announceInterfaces() const noexcept {
  return Scine::Core::DerivedModule::announceInterfaces<OrcaInterfaceModelMap>();
}

std::vector<std::string> OrcaModule::announceModels(const std::string& interface) const noexcept {
  auto models = Scine::Core::DerivedModule::announceModels<OrcaInterfaceModelMap>(interface);

  auto orcaModel = ExternalQC::OrcaCalculator::model;
  if (interface == Scine::Core::Calculator::interface && !orcaFound()) {
    models.erase(std::remove(std::begin(models), std::end(models), orcaModel), std::end(models));
  }

  return models;
}

std::shared_ptr<Scine::Core::Module> OrcaModule::make() {
  return std::make_shared<OrcaModule>();
}

} // namespace Utils
} // namespace Scine
