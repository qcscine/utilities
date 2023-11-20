/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "LennardJonesModule.h"
#include <Core/DerivedModule.h>
#include <Utils/LennardJonesCalculator/LennardJonesCalculator.h>

namespace Scine {
namespace Utils {

using LennardJonesInterfaceModelMap =
    boost::mpl::map<boost::mpl::pair<Scine::Core::Calculator, boost::mpl::vector<LennardJonesCalculator>>>;

std::string LennardJonesModule::name() const noexcept {
  return "Lennardjones";
}

boost::any LennardJonesModule::get(const std::string& interface, const std::string& model) const {
  boost::any resolved = Scine::Core::DerivedModule::resolve<LennardJonesInterfaceModelMap>(interface, model);
  // Throw an exception if we could not match an interface or model
  if (resolved.empty()) {
    throw Scine::Core::ClassNotImplementedError();
  }
  return resolved;
}

bool LennardJonesModule::has(const std::string& interface, const std::string& model) const noexcept {
  if (interface == Scine::Core::Calculator::interface) {
    if (model == LennardJonesCalculator::model) {
      return true;
    }
  }

  return Scine::Core::DerivedModule::has<LennardJonesInterfaceModelMap>(interface, model);
}

std::vector<std::string> LennardJonesModule::announceInterfaces() const noexcept {
  return Scine::Core::DerivedModule::announceInterfaces<LennardJonesInterfaceModelMap>();
}

std::vector<std::string> LennardJonesModule::announceModels(const std::string& interface) const noexcept {
  return Scine::Core::DerivedModule::announceModels<LennardJonesInterfaceModelMap>(interface);
}

std::shared_ptr<Scine::Core::Module> LennardJonesModule::make() {
  return std::make_shared<LennardJonesModule>();
}

} // namespace Utils
} // namespace Scine
