/**
 * @file
 * @brief Default module that is always loaded by Core
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "MrccModule.h"
#include <Core/DerivedModule.h>
#include <Utils/ExternalQC/MRCC/MrccCCCalculator.h>
#include <Utils/ExternalQC/MRCC/MrccDFTCalculator.h>
#include <Utils/ExternalQC/MRCC/MrccHFCalculator.h>
#include <Utils/ExternalQC/MRCC/MrccMP2Calculator.h>

namespace Scine {
namespace Utils {

using MrccInterfaceModelMap =
    boost::mpl::map<boost::mpl::pair<Scine::Core::Calculator, boost::mpl::vector<ExternalQC::MrccHFCalculator, ExternalQC::MrccDFTCalculator,
                                                                                 ExternalQC::MrccMP2Calculator, ExternalQC::MrccCCCalculator>>>;

std::string MrccModule::name() const noexcept {
  return "MRCC";
}

bool MrccModule::mrccFound() {
  return std::getenv(ExternalQC::MrccCalculator::binaryEnvVariable) != nullptr;
}

boost::any MrccModule::get(const std::string& interface, const std::string& model) const {
  boost::any resolved = Scine::Core::DerivedModule::resolve<MrccInterfaceModelMap>(interface, model);
  // Throw an exception if we could not match an interface or model
  if (resolved.empty()) {
    throw Scine::Core::ClassNotImplementedError();
  }
  return resolved;
}

bool MrccModule::has(const std::string& interface, const std::string& model) const noexcept {
  return Scine::Core::DerivedModule::has<MrccInterfaceModelMap>(interface, model) && mrccFound();
}

std::vector<std::string> MrccModule::announceInterfaces() const noexcept {
  return Scine::Core::DerivedModule::announceInterfaces<MrccInterfaceModelMap>();
}

std::vector<std::string> MrccModule::announceModels(const std::string& interface) const noexcept {
  auto models = Scine::Core::DerivedModule::announceModels<MrccInterfaceModelMap>(interface);

  if (interface == Scine::Core::Calculator::interface && !mrccFound()) {
    auto mrccModels = {ExternalQC::MrccHFCalculator::model, ExternalQC::MrccDFTCalculator::model,
                       ExternalQC::MrccMP2Calculator::model, ExternalQC::MrccCCCalculator::model};
    for (const auto& mrccModel : mrccModels) {
      models.erase(std::remove(std::begin(models), std::end(models), mrccModel), std::end(models));
    }
  }
  return models;
}

std::shared_ptr<Scine::Core::Module> MrccModule::make() {
  return std::make_shared<MrccModule>();
}

} // namespace Utils
} // namespace Scine
