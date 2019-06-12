/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_FOCKDIISMODIFIER_H
#define UTILS_FOCKDIISMODIFIER_H

#include "FockDiis.h"
#include <Utils/MethodEssentials/Methods/SCFModifier.h>

namespace Scine {
namespace Utils {

/*!
 * SCFModifier for the DIIS algorithm.
 * This class does not perform directly the DIIS algorithm, it is like an
 * interface between the SCF method and the class for the actual DIIS algorithm.
 */
class FockDiisModifier : public SCFModifier {
 public:
  void onOverlapCalculated() override;
  void onFockCalculated() override;
  void setOrthogonal(bool o) {
    mixer_.setOrthogonal(o);
  }

  void initialize() override;
  void setSpaceSize(int n);

 private:
  bool sameNumberOfElectronsInMethodAndInDensityMatrix();
  FockDiis mixer_;

  bool initialized{false};
};

} // namespace Utils
} // namespace Scine
#endif // FOCKDIISMODIFIER_H
