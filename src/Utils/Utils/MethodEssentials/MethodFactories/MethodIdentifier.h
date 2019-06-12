/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_METHODIDENTIFIER_H
#define UTILS_METHODIDENTIFIER_H

#include <string>

namespace Scine {
namespace Utils {
/*!
 * Identifier for some method type.
 * Fulfills the same role as an enum that can be changed at runtime.
 */
class MethodIdentifier {
 public:
  explicit MethodIdentifier(std::string id = "") : id_(std::move(id)) {
  }
  bool valid() const {
    return !id_.empty();
  }
  bool operator==(const MethodIdentifier& rhs) const {
    return id_ == rhs.id_;
  }
  bool operator!=(const MethodIdentifier& rhs) const {
    return !operator==(rhs);
  }
  bool operator<(const MethodIdentifier& rhs) const {
    return id_ < rhs.id_;
  }

 private:
  std::string id_;
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_METHODIDENTIFIER_H