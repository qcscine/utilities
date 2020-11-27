/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_EDIISMODIFIER_H
#define UTILS_EDIISMODIFIER_H

#include "Ediis.h"
#include <Utils/Scf/MethodInterfaces/ScfModifier.h>

namespace Scine {
namespace Utils {

/*!
 * ScfModifier for the EDIIS algorithm.
 * This class does not perform directly the EDIIS algorithm, it is like an
 * interface between the SCF method and the class for the actual EDIIS algorithm.
 */
class EdiisModifier : public ScfModifier {
 public:
  EdiisModifier() = default;

  void onOverlapCalculated() override;
  void onFockCalculated() override;

  void setSpaceSize(int n);

 private:
  bool sameNumberOfElectronsInMethodAndInDensityMatrix();

  Ediis mixer_;
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_EDIISMODIFIER_H
