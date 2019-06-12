/**
 * @file
 * @brief Default module that is always loaded by Core
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef INCLUDE_SCINE_OS_UTILS_MODULE_H
#define INCLUDE_SCINE_OS_UTILS_MODULE_H

#include "Core/Module.h"
#include "boost/dll/alias.hpp"
#include "boost/hana/define_struct.hpp"
#include <memory>

namespace Scine {
namespace Utils {

/**
 * @brief Module provided by OSUtils. This provides the 'default'
 *   implementations of common interfaces.
 */
class OSUtilsModule : public Core::Module {
 public:
  BOOST_HANA_DEFINE_STRUCT(OSUtilsModule, (std::vector<std::string>, calculator));

  OSUtilsModule() noexcept;

  std::string name() const noexcept final;

  boost::any get(const std::string& interface, const std::string& model) const final;

  bool has(const std::string& interface, const std::string& model) const noexcept final;

  std::vector<std::string> announceInterfaces() const noexcept final;

  std::vector<std::string> announceModels(const std::string& interface) const noexcept final;

  static std::shared_ptr<Module> make();
};

} // namespace Utils
} // namespace Scine

// At global namespace, define the entry point for the module.
BOOST_DLL_ALIAS(Scine::Utils::OSUtilsModule::make, moduleFactory)

#endif