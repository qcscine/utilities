/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Dftd3Parameters.h"
#include "Dftd3ReferencePairs.h"
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/Geometry/ElementTypes.h>
#include <cassert>

namespace Scine {
namespace Utils {
namespace Dftd3 {

// Constructor, which allows you to pass a1, s8 and a2 parameters that belong to the DFT functional used.
Dftd3Parameters::Dftd3Parameters(double a1, double s8, double a2) : a1_(a1), s8_(s8), a2_(a2) {
}

// Wrapper around the covalent radii matrix to provide the correct radius for a given element type.
constexpr double Dftd3Parameters::covalentRadii_[94];
double Dftd3Parameters::getCovalentRadius(ElementType elementType) {
  auto elementTypeIndex = ElementInfo::Z(elementType);
  return covalentRadii_[elementTypeIndex - 1];
}

// Getter for the k1 parameter.
double Dftd3Parameters::getK1() {
  return k1_;
}

// Getter for the k3 parameter.
double Dftd3Parameters::getK3() {
  return k3_;
}

// Wrapper around the D3 reference pair array.
// It returns the 3 x 25 dimensional sub-matrix belonging to the element types of a given atom pair.
constexpr std::array<std::array<std::array<std::array<double, 3>, 25>, 94>, 94> Dftd3ReferencePairs::referencePairs_;
const std::array<std::array<double, 3>, 25>& Dftd3Parameters::getReferencePairs(ElementType elementType1,
                                                                                ElementType elementType2) {
  // Check whether reference pairs are available up to the element types given.
  assert(ElementInfo::Z(elementType1) <= 94 && ElementInfo::Z(elementType2) <= 94);
  return Dftd3ReferencePairs::referencePairs_[ElementInfo::Z(elementType1) - 1][ElementInfo::Z(elementType2) - 1];
}

// Wrapper around the r2r4 parameter matrix.
// Returns the parameter value for a given element type.
constexpr double Dftd3Parameters::r2r4_[94];
double Dftd3Parameters::getR2r4(ElementType elementType) {
  auto elementTypeIndex = ElementInfo::Z(elementType);
  return r2r4_[elementTypeIndex - 1];
}

// Getter for the a1 parameter.
double Dftd3Parameters::getA1() {
  return a1_;
}

// Getter for the s8 parameter.
double Dftd3Parameters::getS8() {
  return s8_;
}

// Getter for the a2 parameter.
double Dftd3Parameters::getA2() {
  return a2_;
}

} // namespace Dftd3
} // namespace Utils
} // namespace Scine
