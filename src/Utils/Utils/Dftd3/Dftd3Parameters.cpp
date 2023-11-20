/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
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

// Wrapper around the r0Zero_ parameter matrix.
// Returns the cut-off radius r0 for a given pair of atoms (only used for zero-damping).
constexpr double Dftd3Parameters::r0Zero_[94][94];
double Dftd3Parameters::getR0Zero(ElementType elementType1, ElementType elementType2) {
  auto elementTypeIndex1 = ElementInfo::Z(elementType1);
  auto elementTypeIndex2 = ElementInfo::Z(elementType2);
  return r0Zero_[elementTypeIndex1 - 1][elementTypeIndex2 - 1];
}

// Getter for the s6 parameter.
double Dftd3Parameters::getS6() {
  return s6_;
}

// Setter for the s6 parameter.
void Dftd3Parameters::setS6(double s6) {
  s6_ = s6;
}

// Getter for the s8 parameter.
double Dftd3Parameters::getS8() {
  return s8_;
}

// Setter for the s8 parameter.
void Dftd3Parameters::setS8(double s8) {
  s8_ = s8;
}

// Getter for the a1 parameter.
double Dftd3Parameters::getA1() {
  return a1_;
}

// Setter for the a1 parameter.
void Dftd3Parameters::setA1(double a1) {
  a1_ = a1;
}

// Getter for the a2 parameter.
double Dftd3Parameters::getA2() {
  return a2_;
}

// Setter for the s6 parameter.
void Dftd3Parameters::setA2(double a2) {
  a2_ = a2;
}

// Getter for the sr parameter.
double Dftd3Parameters::getSr() {
  return sr_;
}

// Setter for the sr parameter.
void Dftd3Parameters::setSr(double sr) {
  sr_ = sr;
}

// Getter for the a parameter.
double Dftd3Parameters::getA() {
  return a_;
}

// Setter for the a parameter.
void Dftd3Parameters::setA(double a) {
  a_ = a;
}

} // namespace Dftd3
} // namespace Utils
} // namespace Scine
