/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_AVAILABLEMETHODS_H
#define UTILS_AVAILABLEMETHODS_H

#include "../MethodFactories/MethodIdentifier.h"
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace Scine {
namespace Utils {

class MethodSpecifier;
class MethodInitializer;
class SinglePointMethod;

/*!
 * This class defines the methods that are available for single-point calculations.
 */
class AvailableMethods {
 public:
  using SpecifierPtr = std::unique_ptr<MethodSpecifier>;
  using SpecifierContainer = std::map<MethodIdentifier, SpecifierPtr>;

  static std::vector<MethodIdentifier> getAvailableMethods();
  static MethodIdentifier addMethod(SpecifierPtr&& specifier);
  static void removeMethod(const MethodIdentifier& specifier);

  static std::string getName(const MethodIdentifier& id);
  static bool isLcaoMethod(const MethodIdentifier& id);
  static bool isScfMethod(const MethodIdentifier& id);
  static bool hasOrthonormalBasis(const MethodIdentifier& id);
  static bool hasNonOrthonormalBasis(const MethodIdentifier& id);
  static unsigned maximalDerivativeOrder(const MethodIdentifier& id);

  /*! Static function to create a RealTimeMethod.
      \param methodID Type of the method to create. */
  static std::unique_ptr<SinglePointMethod> createMethod(const MethodIdentifier& id);
  /*! Given a method, creates the corresponding MethodInitializer.
      Returns nullptr if method is not valid (i.e. of type unknown to this class, or nullptr).*/
  static std::unique_ptr<MethodInitializer> createMethodInitializer(SinglePointMethod* method);

 private:
  static void verifyThatMethodIsRegistered(const MethodIdentifier& id);
  static MethodIdentifier getType(SinglePointMethod* method);
  static SpecifierContainer specifiers_;
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_AVAILABLEMETHODS_H