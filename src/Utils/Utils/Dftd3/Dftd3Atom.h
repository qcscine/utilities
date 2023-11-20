/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SCINE_DFTD3ATOM_H
#define SCINE_DFTD3ATOM_H

#include <Utils/Geometry/Atom.h>

namespace Scine {
namespace Utils {
namespace Dftd3 {
/**
 * @class Dftd3Atom Dftd3Atom.h
 * @brief Describes an atom for a D3 semi-classical dispersion correction calculation.
 */
class Dftd3Atom : public Atom {
 public:
  using Atom::Atom;
  /**
   * @brief Setter for the coordination number (fractional value) for an atom.
   * @param coordinationNumber Coordination number
   */
  void setCoordinationNumber(double coordinationNumber);
  /**
   * @brief Getter for the coordination number (fractional value) for an atom.
   * @return double
   */
  double getCoordinationNumber() const;

  /**
   * @brief Setter for the index of an atom.
   * @param index Index as an integer value.
   */
  void setIndex(int index);
  /**
   * @brief Getter for the index of an atom.
   * @return int
   */
  int getIndex() const;

 private:
  /// @brief coordination number of the atom.
  double coordinationNumber_;
  /// @brief index of the atom.
  int index_;
};

} // namespace Dftd3
} // namespace Utils
} // namespace Scine

#endif // SCINE_DFTD3ATOM_H
