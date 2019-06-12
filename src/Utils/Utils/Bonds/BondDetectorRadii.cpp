/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Bonds/BondDetectorRadii.h"
#include <algorithm>

namespace Scine {
namespace Utils {

BondDetectorRadii::BondDetectorRadii() {
  fillArray();
}

void BondDetectorRadii::fillArray() {
  // initialize with the default value
  std::fill(radii.begin(), radii.end(), double(defaultRadius)); // "double" will not be necessary in future GCC versions

  setRadiusInAngstrom(ElementType::H, 0.23);
  setRadiusInAngstrom(ElementType::He, 0.93);
  setRadiusInAngstrom(ElementType::Li, 0.68);
  setRadiusInAngstrom(ElementType::Be, 0.35);
  setRadiusInAngstrom(ElementType::B, 0.83);
  setRadiusInAngstrom(ElementType::C, 0.68);
  setRadiusInAngstrom(ElementType::N, 0.68);
  setRadiusInAngstrom(ElementType::O, 0.68);
  setRadiusInAngstrom(ElementType::F, 0.64);
  setRadiusInAngstrom(ElementType::Ne, 1.12);
  setRadiusInAngstrom(ElementType::Na, 0.97);
  setRadiusInAngstrom(ElementType::Mg, 1.1);
  setRadiusInAngstrom(ElementType::Al, 1.35);
  setRadiusInAngstrom(ElementType::Si, 1.2);
  setRadiusInAngstrom(ElementType::P, 1.05);
  setRadiusInAngstrom(ElementType::S, 1.02);
  setRadiusInAngstrom(ElementType::Cl, 0.99);
  setRadiusInAngstrom(ElementType::Ar, 1.57);
  setRadiusInAngstrom(ElementType::K, 1.33);
  setRadiusInAngstrom(ElementType::Ca, 0.99);
  setRadiusInAngstrom(ElementType::Sc, 1.44);
  setRadiusInAngstrom(ElementType::Ti, 1.47);
  setRadiusInAngstrom(ElementType::V, 1.33);
  setRadiusInAngstrom(ElementType::Cr, 1.35);
  setRadiusInAngstrom(ElementType::Mn, 1.35);
  setRadiusInAngstrom(ElementType::Fe, 1.34);
  setRadiusInAngstrom(ElementType::Co, 1.33);
  setRadiusInAngstrom(ElementType::Ni, 1.5);
  setRadiusInAngstrom(ElementType::Cu, 1.52);
  setRadiusInAngstrom(ElementType::Zn, 1.45);
  setRadiusInAngstrom(ElementType::Ga, 1.22);
  setRadiusInAngstrom(ElementType::Ge, 1.17);
  setRadiusInAngstrom(ElementType::As, 1.21);
  setRadiusInAngstrom(ElementType::Se, 1.22);
  setRadiusInAngstrom(ElementType::Br, 1.21);
  setRadiusInAngstrom(ElementType::Kr, 1.91);
  setRadiusInAngstrom(ElementType::Rb, 1.47);
  setRadiusInAngstrom(ElementType::Sr, 1.12);
  setRadiusInAngstrom(ElementType::Y, 1.78);
  setRadiusInAngstrom(ElementType::Zr, 1.56);
  setRadiusInAngstrom(ElementType::Nb, 1.48);
  setRadiusInAngstrom(ElementType::Mo, 1.47);
  setRadiusInAngstrom(ElementType::Tc, 1.35);
  setRadiusInAngstrom(ElementType::Ru, 1.4);
  setRadiusInAngstrom(ElementType::Rh, 1.45);
  setRadiusInAngstrom(ElementType::Pd, 1.5);
  setRadiusInAngstrom(ElementType::Ag, 1.59);
  setRadiusInAngstrom(ElementType::Cd, 1.69);
  setRadiusInAngstrom(ElementType::In, 1.63);
  setRadiusInAngstrom(ElementType::Sn, 1.46);
  setRadiusInAngstrom(ElementType::Te, 1.47);
  setRadiusInAngstrom(ElementType::I, 1.4);
  setRadiusInAngstrom(ElementType::Xe, 1.98);
  setRadiusInAngstrom(ElementType::Cs, 1.67);
  setRadiusInAngstrom(ElementType::Ba, 1.34);
  setRadiusInAngstrom(ElementType::La, 1.87);
  setRadiusInAngstrom(ElementType::Ce, 1.83);
  setRadiusInAngstrom(ElementType::Pr, 1.82);
  setRadiusInAngstrom(ElementType::Nd, 1.81);
  setRadiusInAngstrom(ElementType::Pm, 1.8);
  setRadiusInAngstrom(ElementType::Sm, 1.8);
  setRadiusInAngstrom(ElementType::Eu, 1.99);
  setRadiusInAngstrom(ElementType::Gd, 1.79);
  setRadiusInAngstrom(ElementType::Tb, 1.76);
  setRadiusInAngstrom(ElementType::Dy, 1.75);
  setRadiusInAngstrom(ElementType::Ho, 1.74);
  setRadiusInAngstrom(ElementType::Er, 1.73);
  setRadiusInAngstrom(ElementType::Tm, 1.72);
}

void BondDetectorRadii::setRadiusInAngstrom(ElementType e, double rInAngstrom) {
  setRadius(e, toBohr(Angstrom(rInAngstrom)));
}

void BondDetectorRadii::setRadius(ElementType e, double r) {
  auto index = static_cast<int>(e);
  radii[index] = r;
}

} /* namespace Utils */
} /* namespace Scine */
