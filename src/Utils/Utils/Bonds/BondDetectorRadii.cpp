/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Bonds/BondDetectorRadii.h"
#include "Utils/Geometry/ElementInfo.h"
#include <algorithm>

namespace Scine {
namespace Utils {

BondDetectorRadii::BondDetectorRadii() {
  fillArray();
}

void BondDetectorRadii::fillArray() {
  // initialize with the default value
  std::fill(radii.begin(), radii.end(), double(defaultRadius)); // "double" will not be necessary in future GCC versions

  // Fill in values obtained from the Cambridge Structural Database
  // https://www.ccdc.cam.ac.uk/support-and-resources/ccdcresources/Elemental_Radii.xlsx
  setRadiusInAngstrom(ElementType::H, 0.23);
  setRadiusInAngstrom(ElementType::He, 1.5);
  setRadiusInAngstrom(ElementType::Li, 1.28);
  setRadiusInAngstrom(ElementType::Be, 0.96);
  setRadiusInAngstrom(ElementType::B, 0.83);
  setRadiusInAngstrom(ElementType::C, 0.68);
  setRadiusInAngstrom(ElementType::N, 0.68);
  setRadiusInAngstrom(ElementType::O, 0.68);
  setRadiusInAngstrom(ElementType::F, 0.64);
  setRadiusInAngstrom(ElementType::Ne, 1.5);
  setRadiusInAngstrom(ElementType::Na, 1.66);
  setRadiusInAngstrom(ElementType::Mg, 1.41);
  setRadiusInAngstrom(ElementType::Al, 1.21);
  setRadiusInAngstrom(ElementType::Si, 1.2);
  setRadiusInAngstrom(ElementType::P, 1.05);
  setRadiusInAngstrom(ElementType::S, 1.02);
  setRadiusInAngstrom(ElementType::Cl, 0.99);
  setRadiusInAngstrom(ElementType::Ar, 1.51);
  setRadiusInAngstrom(ElementType::K, 2.03);
  setRadiusInAngstrom(ElementType::Ca, 1.76);
  setRadiusInAngstrom(ElementType::Sc, 1.7);
  setRadiusInAngstrom(ElementType::Ti, 1.6);
  setRadiusInAngstrom(ElementType::V, 1.53);
  setRadiusInAngstrom(ElementType::Cr, 1.39);
  setRadiusInAngstrom(ElementType::Mn, 1.61);
  setRadiusInAngstrom(ElementType::Fe, 1.52);
  setRadiusInAngstrom(ElementType::Co, 1.26);
  setRadiusInAngstrom(ElementType::Ni, 1.24);
  setRadiusInAngstrom(ElementType::Cu, 1.32);
  setRadiusInAngstrom(ElementType::Zn, 1.22);
  setRadiusInAngstrom(ElementType::Ga, 1.22);
  setRadiusInAngstrom(ElementType::Ge, 1.17);
  setRadiusInAngstrom(ElementType::As, 1.21);
  setRadiusInAngstrom(ElementType::Se, 1.22);
  setRadiusInAngstrom(ElementType::Br, 1.21);
  setRadiusInAngstrom(ElementType::Kr, 1.5);
  setRadiusInAngstrom(ElementType::Rb, 2.2);
  setRadiusInAngstrom(ElementType::Sr, 1.95);
  setRadiusInAngstrom(ElementType::Y, 1.9);
  setRadiusInAngstrom(ElementType::Zr, 1.75);
  setRadiusInAngstrom(ElementType::Nb, 1.64);
  setRadiusInAngstrom(ElementType::Mo, 1.54);
  setRadiusInAngstrom(ElementType::Tc, 1.47);
  setRadiusInAngstrom(ElementType::Ru, 1.46);
  setRadiusInAngstrom(ElementType::Rh, 1.42);
  setRadiusInAngstrom(ElementType::Pd, 1.39);
  setRadiusInAngstrom(ElementType::Ag, 1.45);
  setRadiusInAngstrom(ElementType::Cd, 1.54);
  setRadiusInAngstrom(ElementType::In, 1.42);
  setRadiusInAngstrom(ElementType::Sn, 1.39);
  setRadiusInAngstrom(ElementType::Sb, 1.39);
  setRadiusInAngstrom(ElementType::Te, 1.47);
  setRadiusInAngstrom(ElementType::I, 1.4);
  setRadiusInAngstrom(ElementType::Xe, 1.5);
  setRadiusInAngstrom(ElementType::Cs, 2.44);
  setRadiusInAngstrom(ElementType::Ba, 2.15);
  setRadiusInAngstrom(ElementType::La, 2.07);
  setRadiusInAngstrom(ElementType::Ce, 2.04);
  setRadiusInAngstrom(ElementType::Pr, 2.03);
  setRadiusInAngstrom(ElementType::Nd, 2.01);
  setRadiusInAngstrom(ElementType::Pm, 1.99);
  setRadiusInAngstrom(ElementType::Sm, 1.98);
  setRadiusInAngstrom(ElementType::Eu, 1.98);
  setRadiusInAngstrom(ElementType::Gd, 1.96);
  setRadiusInAngstrom(ElementType::Tb, 1.94);
  setRadiusInAngstrom(ElementType::Dy, 1.92);
  setRadiusInAngstrom(ElementType::Ho, 1.92);
  setRadiusInAngstrom(ElementType::Er, 1.89);
  setRadiusInAngstrom(ElementType::Tm, 1.9);
  setRadiusInAngstrom(ElementType::Yb, 1.87);
  setRadiusInAngstrom(ElementType::Lu, 1.87);
  setRadiusInAngstrom(ElementType::Hf, 1.75);
  setRadiusInAngstrom(ElementType::Ta, 1.7);
  setRadiusInAngstrom(ElementType::W, 1.62);
  setRadiusInAngstrom(ElementType::Re, 1.51);
  setRadiusInAngstrom(ElementType::Os, 1.44);
  setRadiusInAngstrom(ElementType::Ir, 1.41);
  setRadiusInAngstrom(ElementType::Pt, 1.36);
  setRadiusInAngstrom(ElementType::Au, 1.36);
  setRadiusInAngstrom(ElementType::Hg, 1.32);
  setRadiusInAngstrom(ElementType::Tl, 1.45);
  setRadiusInAngstrom(ElementType::Pb, 1.46);
  setRadiusInAngstrom(ElementType::Bi, 1.48);
  setRadiusInAngstrom(ElementType::Po, 1.4);
  setRadiusInAngstrom(ElementType::At, 1.21);
  setRadiusInAngstrom(ElementType::Rn, 1.5);
  setRadiusInAngstrom(ElementType::Fr, 2.6);
  setRadiusInAngstrom(ElementType::Ra, 2.21);
  setRadiusInAngstrom(ElementType::Ac, 2.15);
  setRadiusInAngstrom(ElementType::Th, 2.06);
  setRadiusInAngstrom(ElementType::Pa, 2.0);
  setRadiusInAngstrom(ElementType::U, 1.96);
  setRadiusInAngstrom(ElementType::Np, 1.9);
  setRadiusInAngstrom(ElementType::Pu, 1.87);
  setRadiusInAngstrom(ElementType::Am, 1.8);
  setRadiusInAngstrom(ElementType::Cm, 1.69);
  setRadiusInAngstrom(ElementType::Bk, 1.54);
  setRadiusInAngstrom(ElementType::Cf, 1.83);
  setRadiusInAngstrom(ElementType::Es, 1.5);
  setRadiusInAngstrom(ElementType::Fm, 1.5);
  setRadiusInAngstrom(ElementType::Md, 1.5);
  setRadiusInAngstrom(ElementType::No, 1.5);
  setRadiusInAngstrom(ElementType::Rf, 1.5);
  setRadiusInAngstrom(ElementType::Db, 1.5);
  setRadiusInAngstrom(ElementType::Sg, 1.5);
  setRadiusInAngstrom(ElementType::Bh, 1.5);
  setRadiusInAngstrom(ElementType::Hs, 1.5);
  setRadiusInAngstrom(ElementType::Mt, 1.5);
  setRadiusInAngstrom(ElementType::Ds, 1.5);
}

double BondDetectorRadii::getRadius(ElementType e) const {
  auto index = ElementInfo::Z(e);
  return radii[index];
}

void BondDetectorRadii::setRadiusInAngstrom(ElementType e, double rInAngstrom) {
  setRadius(e, toBohr(Angstrom(rInAngstrom)));
}

void BondDetectorRadii::setRadius(ElementType e, double r) {
  auto index = ElementInfo::Z(e);
  radii[index] = r;
}

} /* namespace Utils */
} /* namespace Scine */
