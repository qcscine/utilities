/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_CHARGEMODEL5_H
#define UTILS_CHARGEMODEL5_H

#include <Utils/Typenames.h>
#include <vector>

namespace Scine {
namespace Utils {
class AtomCollection;

/**
 * @class ChargeModel5 ChargeModel5.h
 * @brief This class provides the functionality to obtain Charge Model 5 charges from Hirshfeld charges.
 *
 *        Implemented based on this paper: J. Chem. Theory Comput. 2012, 8, 2, 527-541.
 */
class ChargeModel5 {
 public:
  /**
   * @brief Takes a set of Hirshfeld charges and a molecular structure as an input and returns the CM5 charges.
   * @param hirshfeldCharges The Hirshfeld charges.
   * @param atoms The molecular structure.
   * @return The CM5 charges.
   */
  static std::vector<double> calculateCm5Charges(const std::vector<double>& hirshfeldCharges, const AtomCollection& atoms);

 private:
  // Getter for the pairwise CM5 parameters
  static double getPairwiseParameter(const ElementType& e1, const ElementType& e2);
  // The global alpha parameter
  static constexpr double globalAlphaParameter_ = 2.474;
  // The atomwise parameters for the CM5 charges
  static constexpr double atomwiseParameters_[118] = {
      0.0056,  -0.1543, 0.0,     0.0333,  -0.103,  -0.0446, -0.1072, -0.0802, -0.0629, -0.1088, 0.0184,  0.0,
      -0.0726, -0.079,  -0.0756, -0.0565, -0.0444, -0.0767, 0.013,   0.0,     0.0,     0.0,     0.0,     0.0,
      0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     -0.0512, -0.0557, -0.0533, -0.0399, -0.0313, -0.0541,
      0.0092,  0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
      -0.0361, -0.0393, -0.0376, -0.0281, -0.022,  -0.0381, 0.0065,  0.0,     0.0,     0.0,     0.0,     0.0,
      0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
      0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     -0.0255, -0.0277, -0.0265, -0.0198,
      -0.0155, -0.0269, 0.0046,  0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
      0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
      0.0,     0.0,     0.0,     0.0,     -0.0179, -0.0195, -0.0187, -0.014,  -0.011,  -0.0189};
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_CHARGEMODEL5_H
