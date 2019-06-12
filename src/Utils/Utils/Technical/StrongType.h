/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_STRONGTYPE_H_
#define UTILS_STRONGTYPE_H_

#include <utility>

namespace Scine {
namespace Utils {

/**
 * @class StrongType StrongType.h
 * From: https://www.fluentcpp.com/2016/12/08/strong-types-for-strong-interfaces/
 * @tparam T The underlying type.
 * @tparam Parameter A "phantom type" serves the purpose of specializing the type.
 */
template<typename T, typename Parameter>
class StrongType {
 public:
  /**
   * @brief Copy Constructor.
   */
  constexpr explicit StrongType(T const& value) : value_(value) {
  }
  /**
   * @brief Move Constructor.
   */
  constexpr explicit StrongType(T&& value) : value_(std::move(value)) {
  }
  /**
   * @brief Getter for the underlying base class.
   * @return The value in its base from.
   */
  constexpr const T& get() const {
    return value_;
  }

 private:
  T value_;
};

} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_STRONGTYPE_H_
