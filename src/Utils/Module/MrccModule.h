/**
 * @file
 * @brief Default module that is always loaded by Core
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_MRCCMODULE_H
#define UTILS_MRCCMODULE_H

#include "Core/Module.h"

namespace Scine {
namespace Utils {
/**
 * @class
 * @brief MRCC module provided by OSUtils.
 */
class MrccModule : public Core::Module {
 public:
  std::string name() const noexcept final;

  boost::any get(const std::string& interface, const std::string& model) const final;

  bool has(const std::string& interface, const std::string& model) const noexcept final;

  std::vector<std::string> announceInterfaces() const noexcept final;

  std::vector<std::string> announceModels(const std::string& interface) const noexcept final;

  static std::shared_ptr<Core::Module> make();

  static bool mrccFound();
};

// Shared library entry point creating pointers to all contained modules
std::vector<std::shared_ptr<Core::Module>> moduleFactory();

} // namespace Utils
} // namespace Scine

#endif // UTILS_MRCCMODULE_H
