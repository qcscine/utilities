/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_FOCKDIISMODIFIER_H
#define UTILS_FOCKDIISMODIFIER_H

#include "FockDiis.h"
#include <Utils/Scf/MethodInterfaces/ScfModifier.h>

namespace Scine {
namespace Utils {

/*!
 * ScfModifier for the DIIS algorithm.
 * This class does not perform directly the DIIS algorithm, it is like an
 * interface between the Scf method and the class for the actual DIIS algorithm.
 */
class FockDiisModifier : public ScfModifier {
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
