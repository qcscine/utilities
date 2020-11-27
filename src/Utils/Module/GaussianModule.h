/**
 * @file
 * @brief Default module that is always loaded by Core
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SCINE_GAUSSIAN_MODULE_H
#define SCINE_GAUSSIAN_MODULE_H

#include "Core/Module.h"
#include "boost/dll/alias.hpp"
#include <memory>

namespace Scine {
namespace Utils {

/**
 * @brief Gaussian module provided by OSUtils.
 */
class GaussianModule : public Core::Module {
 public:
  static bool gaussianFound();

  std::string name() const noexcept final;

  boost::any get(const std::string& interface, const std::string& model) const final;

  bool has(const std::string& interface, const std::string& model) const noexcept final;

  std::vector<std::string> announceInterfaces() const noexcept final;

  std::vector<std::string> announceModels(const std::string& interface) const noexcept final;

  static std::shared_ptr<Core::Module> make();
};

} // namespace Utils
} // namespace Scine

#endif // SCINE_GAUSSIAN_MODULE_H
