/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "AvailableMethods.h"
#include "MethodInitializer.h"
#include "MethodSpecifier.h"
#include <Utils/MethodEssentials/Methods/SinglePointMethod.h>

namespace Scine {
namespace Utils {

AvailableMethods::SpecifierContainer AvailableMethods::specifiers_;

std::vector<MethodIdentifier> AvailableMethods::getAvailableMethods() {
  std::vector<MethodIdentifier> methods;
  for (const auto& s : specifiers_)
    methods.push_back(s.first);
  return methods;
}

MethodIdentifier AvailableMethods::addMethod(AvailableMethods::SpecifierPtr&& specifier) {
  MethodIdentifier identifier(specifier->getName());
  specifiers_.emplace(identifier, std::move(specifier));
  return identifier;
}

void AvailableMethods::removeMethod(const MethodIdentifier& specifier) {
  specifiers_.erase(specifier);
}

std::string AvailableMethods::getName(const MethodIdentifier& id) {
  verifyThatMethodIsRegistered(id);
  return specifiers_[id]->getName();
}

bool AvailableMethods::isLcaoMethod(const MethodIdentifier& id) {
  verifyThatMethodIsRegistered(id);
  return specifiers_[id]->isLcaoMethod();
}

bool AvailableMethods::isScfMethod(const MethodIdentifier& id) {
  verifyThatMethodIsRegistered(id);
  return specifiers_[id]->isScfMethod();
}

bool AvailableMethods::hasOrthonormalBasis(const MethodIdentifier& id) {
  verifyThatMethodIsRegistered(id);
  return specifiers_[id]->hasOrthonormalBasis();
}

bool AvailableMethods::hasNonOrthonormalBasis(const MethodIdentifier& id) {
  return !hasOrthonormalBasis(id);
}

unsigned AvailableMethods::maximalDerivativeOrder(const MethodIdentifier& id) {
  verifyThatMethodIsRegistered(id);
  return specifiers_[id]->maximalDerivativeOrder();
}

void AvailableMethods::verifyThatMethodIsRegistered(const MethodIdentifier& id) {
  if (specifiers_.find(id) == specifiers_.end())
    throw std::runtime_error("Requested method is not registered.");
}

std::unique_ptr<SinglePointMethod> AvailableMethods::createMethod(const MethodIdentifier& id) {
  verifyThatMethodIsRegistered(id);
  return specifiers_[id]->createMethod();
}

std::unique_ptr<MethodInitializer> AvailableMethods::createMethodInitializer(SinglePointMethod* method) {
  if (method == nullptr)
    return nullptr;
  auto methodID = getType(method);
  verifyThatMethodIsRegistered(methodID);
  return specifiers_[methodID]->createInitializer(method);
}

MethodIdentifier AvailableMethods::getType(SinglePointMethod* method) {
  for (const auto& specifier : specifiers_) {
    if (specifier.second->compatibleType(method)) {
      return specifier.first;
    }
  }
  throw std::runtime_error("No MethodIdentifier is compatible with the given pointer");
}
} // namespace Utils
} // namespace Scine
