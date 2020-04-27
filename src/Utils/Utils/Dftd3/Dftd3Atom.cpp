/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Dftd3Atom.h"

namespace Scine {
namespace Utils {
namespace Dftd3 {

// Setter for the coordination number of an atom.
void Dftd3Atom::setCoordinationNumber(double coordinationNumber) {
  coordinationNumber_ = coordinationNumber;
}

// Getter for the coordination number of an atom.
double Dftd3Atom::getCoordinationNumber() const {
  return coordinationNumber_;
}

// Setter for the index of an atom.
void Dftd3Atom::setIndex(int index) {
  index_ = index;
}

// Getter for the index of an atom.
int Dftd3Atom::getIndex() const {
  return index_;
}

} // namespace Dftd3
} // namespace Utils
} // namespace Scine