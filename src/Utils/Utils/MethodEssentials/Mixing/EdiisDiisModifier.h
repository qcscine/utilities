/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_EDIISDIISMODIFIER_H
#define UTILS_EDIISDIISMODIFIER_H

#include "Ediis.h"
#include "FockDiis.h"
#include <Utils/MethodEssentials/Methods/SCFModifier.h>

namespace Scine {
namespace Utils {

/*!
 * Class to use the EDIIS + DIIS combination as in Gaussian.
 * Garza, Scuseria, J Chem Phys 137 (2012) 054110:
 * "In the default EDIIS + DIIS procedure in Gaussian, EDIIS is used when the largest DIIS error (errMax) is greater
 * than 10^-1 a.u. but DIIS is employed when this error goes below 10^-4 a.u. In between these values, the EDIIS and
 * DIIS coefficients are weighted such that c = 10(errorMax) C(EDIIS) + (1-10(errorMax))c(DIIS); however, if the error
 * of the last cycle is 10% greater than the minimum error, pure EDIIS is used."
 */
class EdiisDiisModifier : public SCFModifier {
 public:
  EdiisDiisModifier();
  void onOverlapCalculated() override;
  void onFockCalculated() override;
  void setSpaceSize(unsigned n);
  void initialize() override;

 private:
  bool sameNumberOfElectronsInMethodAndInDensityMatrix();
  void setOrthogonal(bool o);

  FockDiis diis_;
  Ediis ediis_;

  bool initialized;
  SpinAdaptedMatrix getCombinedFockMatrix();
  SpinAdaptedMatrix mixedFockMatrix(double errMax);
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_EDIISDIISMODIFIER_H
