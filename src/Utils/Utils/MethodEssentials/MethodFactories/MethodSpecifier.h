/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_METHODSPECIFIER_H
#define UTILS_METHODSPECIFIER_H

#include <memory>
#include <string>

namespace Scine {
namespace Utils {

class SinglePointMethod;
class MethodInitializer;

// Determines whether m is of type T*
template<typename T>
inline bool checkMethodType(SinglePointMethod* m) {
  return dynamic_cast<T*>(m) != nullptr;
}

/*!
 * This class defines the properties of a method and implements a factory to create instances of it.
 */
class MethodSpecifier {
 public:
  virtual ~MethodSpecifier() = default;

  virtual std::string getName() const = 0;
  virtual bool isLcaoMethod() const = 0;
  virtual bool isScfMethod() const = 0;
  virtual bool hasOrthonormalBasis() const = 0;
  virtual unsigned maximalDerivativeOrder() const = 0;
  virtual bool compatibleType(SinglePointMethod* method) const = 0;
  virtual std::unique_ptr<SinglePointMethod> createMethod() const = 0;
  virtual std::unique_ptr<MethodInitializer> createInitializer(SinglePointMethod* method) const = 0;
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_METHODSPECIFIER_H