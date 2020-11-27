/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_SCOPEDLOCALE_H_
#define UTILS_SCOPEDLOCALE_H_

#include <string>

namespace Scine {
namespace Utils {

/**
 * @class ScopedLocale ScopedLocale.h
 * @brief Introduces a scope to change the locale, the original locale will be set back upon destruction.
 */
class ScopedLocale {
 public:
  /// @brief Explicit constructor to avoid unintentional instances.
  explicit ScopedLocale(const std::string& targetLocale);
  /// @brief Deleted copy constructor to avoid duplicate instances.
  ScopedLocale(const ScopedLocale&) = delete;
  /// @brief Move constructor.
  ScopedLocale(ScopedLocale&& /*other*/) noexcept;
  /// @brief Custom Destructor.
  ~ScopedLocale();
  ScopedLocale& operator=(const ScopedLocale&) = delete;
  ScopedLocale& operator=(ScopedLocale&& /*other*/) noexcept;

  /**
   * @brief Scope for the standard C locale.
   * @return ScopedLocale
   */
  static ScopedLocale cLocale();

 private:
  std::string originalLocale_;
  bool mustBeSetBackInDestructor_ = true;
};

} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_SCOPEDLOCALE_H_
