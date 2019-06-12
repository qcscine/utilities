/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_ATOMICGTOS_H
#define UTILS_ATOMICGTOS_H

#include <Utils/MethodEssentials/util/GTOExpansion.h>
#include <utility>

namespace Scine {
namespace Utils {

/*!
 * Class containing the GTO expansions for
 * one atom.
 */

class AtomicGTOs {
 public:
  bool hasS() const {
    return validS_;
  }

  bool hasP() const {
    return validP_;
  }

  bool hasD() const {
    return validD_;
  }

  void setS(GTOExpansion g) {
    s_ = std::move(g);
    validS_ = true;
  }

  void setP(GTOExpansion g) {
    p_ = std::move(g);
    validP_ = true;
  }

  void setD(GTOExpansion g) {
    d_ = std::move(g);
    validD_ = true;
  }

  const GTOExpansion& s() const {
    return s_;
  }

  const GTOExpansion& p() const {
    return p_;
  }

  const GTOExpansion& d() const {
    return d_;
  }

 private:
  bool validS_{false}, validP_{false}, validD_{false};
  GTOExpansion s_, p_, d_;
};
} // namespace Utils
} // namespace Scine
#endif // UTILS_ATOMICGTOS_H
