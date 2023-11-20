/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "TurbomoleModule.h"
#include <Core/DerivedModule.h>
#include <Utils/ExternalQC/Turbomole/TurbomoleCalculator.h>

namespace Scine {
namespace Utils {

using TurbomoleInterfaceModelMap =
    boost::mpl::map<boost::mpl::pair<Scine::Core::Calculator, boost::mpl::vector<ExternalQC::TurbomoleCalculator>>>;

std::string TurbomoleModule::name() const noexcept {
  return "Turbomole";
}

bool TurbomoleModule::turbomoleFound() {
  return std::getenv("TURBODIR") != nullptr;
}

boost::any TurbomoleModule::get(const std::string& interface, const std::string& model) const {
  if (interface == Scine::Core::Calculator::interface && model == ExternalQC::TurbomoleCalculator::model && !turbomoleFound()) {
    throw Scine::Core::ClassNotImplementedError();
  }

  boost::any resolved = Scine::Core::DerivedModule::resolve<TurbomoleInterfaceModelMap>(interface, model);
  // Throw an exception if we could not match an interface or model
  if (resolved.empty()) {
    throw Scine::Core::ClassNotImplementedError();
  }
  return resolved;
}

bool TurbomoleModule::has(const std::string& interface, const std::string& model) const noexcept {
  if (interface == Scine::Core::Calculator::interface) {
    if (model == ExternalQC::TurbomoleCalculator::model)
      return turbomoleFound();
  }

  return Scine::Core::DerivedModule::has<TurbomoleInterfaceModelMap>(interface, model);
}

std::vector<std::string> TurbomoleModule::announceInterfaces() const noexcept {
  return Scine::Core::DerivedModule::announceInterfaces<TurbomoleInterfaceModelMap>();
}

std::vector<std::string> TurbomoleModule::announceModels(const std::string& interface) const noexcept {
  auto models = Scine::Core::DerivedModule::announceModels<TurbomoleInterfaceModelMap>(interface);

  auto turbomoleModel = ExternalQC::TurbomoleCalculator::model;
  if (interface == Scine::Core::Calculator::interface && !turbomoleFound()) {
    models.erase(std::remove(std::begin(models), std::end(models), turbomoleModel), std::end(models));
  }

  return models;
}

std::shared_ptr<Scine::Core::Module> TurbomoleModule::make() {
  return std::make_shared<TurbomoleModule>();
}

} // namespace Utils
} // namespace Scine
