/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "ChargeModel5.h"
#include <Utils/Constants.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Geometry/ElementInfo.h>

namespace Scine {
namespace Utils {

std::vector<double> ChargeModel5::calculateCm5Charges(const std::vector<double>& hirshfeldCharges, const AtomCollection& atoms) {
  if (hirshfeldCharges.size() != atoms.size())
    throw std::runtime_error("The number of atoms is not the same as the size of the Hirshfeld charges vector.");

  std::vector<double> cm5Charges;
  for (int i = 0; i < atoms.size(); ++i) {
    double charge = hirshfeldCharges[i];
    double covRad = ElementInfo::covalentRadius(atoms.getElement(i));
    for (int j = 0; j < atoms.size(); ++j) {
      if (i == j)
        continue;
      double distance = (atoms.getPosition(i) - atoms.getPosition(j)).norm();
      double paulingBondOrder = std::exp(-globalAlphaParameter_ * Constants::angstrom_per_bohr *
                                         (distance - covRad - ElementInfo::covalentRadius(atoms.getElement(j))));
      charge += paulingBondOrder * getPairwiseParameter(atoms.getElement(i), atoms.getElement(j));
    }
    cm5Charges.push_back(charge);
  }
  return cm5Charges;
}

constexpr double ChargeModel5::atomwiseParameters_[118];
double ChargeModel5::getPairwiseParameter(const ElementType& e1, const ElementType& e2) {
  auto z1 = ElementInfo::Z(e1);
  auto z2 = ElementInfo::Z(e2);

  if (z1 == z2) // implicitly already covered by last case
    return 0.0;
  else if (z1 == 1 && z2 == 6)
    return 0.0502;
  else if (z1 == 6 && z2 == 1)
    return -0.0502;
  else if (z1 == 1 && z2 == 7)
    return 0.1747;
  else if (z1 == 7 && z2 == 1)
    return -0.1747;
  else if (z1 == 1 && z2 == 8)
    return 0.1671;
  else if (z1 == 8 && z2 == 1)
    return -0.1671;
  else if (z1 == 6 && z2 == 7)
    return 0.0556;
  else if (z1 == 7 && z2 == 6)
    return -0.0556;
  else if (z1 == 6 && z2 == 8)
    return 0.0234;
  else if (z1 == 8 && z2 == 6)
    return -0.0234;
  else if (z1 == 7 && z2 == 8)
    return -0.0346;
  else if (z1 == 8 && z2 == 7)
    return 0.0346;
  else
    return atomwiseParameters_[z1 - 1] - atomwiseParameters_[z2 - 1];
}

} // namespace Utils
} // namespace Scine
