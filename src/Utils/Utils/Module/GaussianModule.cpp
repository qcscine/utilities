/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Core/DerivedModule.h>
#include <Utils/ExternalQC/Gaussian/GaussianCalculator.h>
#include <Utils/Module/GaussianModule.h>

namespace Scine {
namespace Utils {

using InterfaceModelMap =
    boost::mpl::map<boost::mpl::pair<Scine::Core::Calculator, boost::mpl::vector<ExternalQC::GaussianCalculator>>>;

std::string GaussianModule::name() const noexcept {
  return "Gaussian";
};

bool GaussianModule::gaussianFound() {
  return std::getenv("GAUSSIAN_BINARY_PATH") != nullptr;
}

boost::any GaussianModule::get(const std::string& interface, const std::string& model) const {
  if (interface == Scine::Core::Calculator::interface && model == ExternalQC::GaussianCalculator::model && !gaussianFound()) {
    throw Scine::Core::ClassNotImplementedError();
  }

  boost::any resolved = Scine::Core::DerivedModule::resolve<InterfaceModelMap>(interface, model);
  // Throw an exception if we could not match an interface or model
  if (resolved.empty()) {
    throw Scine::Core::ClassNotImplementedError();
  }
  return resolved;
}

bool GaussianModule::has(const std::string& interface, const std::string& model) const noexcept {
  if (interface == Scine::Core::Calculator::interface) {
    if (model == ExternalQC::GaussianCalculator::model) {
      return gaussianFound();
    }
  }

  return Scine::Core::DerivedModule::has<InterfaceModelMap>(interface, model);
}

std::vector<std::string> GaussianModule::announceInterfaces() const noexcept {
  return Scine::Core::DerivedModule::announceInterfaces<InterfaceModelMap>();
}

std::vector<std::string> GaussianModule::announceModels(const std::string& interface) const noexcept {
  auto models = Scine::Core::DerivedModule::announceModels<InterfaceModelMap>(interface);

  auto gaussianModel = ExternalQC::GaussianCalculator::model;
  if (interface == Scine::Core::Calculator::interface && !gaussianFound()) {
    models.erase(std::remove(std::begin(models), std::end(models), gaussianModel), std::end(models));
  }

  return models;
}

std::shared_ptr<Scine::Core::Module> GaussianModule::make() {
  return std::make_shared<GaussianModule>();
}

} // namespace Utils
} // namespace Scine
