/**
 * @file
 * @brief Default module that is always loaded by Core
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SCINE_TURBOMOLE_MODULE_H
#define SCINE_TURBOMOLE_MODULE_H

#include "Core/Module.h"
#include "boost/dll/alias.hpp"
#include <memory>

namespace Scine {
namespace Utils {

/**
 * @brief Turbomole Module provided by OSUtils.
 */
class TurbomoleModule : public Core::Module {
 public:
  static bool turbomoleFound();

  std::string name() const noexcept final;

  boost::any get(const std::string& interface, const std::string& model) const final;

  bool has(const std::string& interface, const std::string& model) const noexcept final;

  std::vector<std::string> announceInterfaces() const noexcept final;

  std::vector<std::string> announceModels(const std::string& interface) const noexcept final;

  static std::shared_ptr<Core::Module> make();
};

// Shared library entry point creating pointers to all contained modules
std::vector<std::shared_ptr<Core::Module>> moduleFactory();

} // namespace Utils
} // namespace Scine

#endif // SCINE_TURBOMOLE_MODULE_H