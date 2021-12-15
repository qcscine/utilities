/**
 * @file ElementData.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Geometry/ElementData.h"
#include "Utils/Geometry/ElementInfo.h"

namespace Scine {
namespace Utils {
namespace Constants {

int ElementDataSingleton::ElementData::valElectrons() const {
  // Possibility to specify exceptions
  if (d_valElectrons > -1) {
    return d_valElectrons;
  }

  const bool isMainGroup = (d_Z <= 20 || (31 <= d_Z && d_Z <= 38) || (49 <= d_Z && d_Z <= 56) ||
                            (81 <= d_Z && d_Z <= 88) || (113 <= d_Z && d_Z <= 118));

  if (isMainGroup) {
    return d_sElectrons + d_pElectrons;
  }

  const bool isFGroup = ((57 <= d_Z && d_Z <= 70) || (89 <= d_Z && d_Z <= 102));

  if (isFGroup) {
    return d_sElectrons + d_fElectrons;
  }

  return d_sElectrons + d_dElectrons;
}

ElementDataSingleton::ElementData ElementDataSingleton::operator[](ElementType type) const {
  return lookup(type);
}

ElementDataSingleton::ElementData ElementDataSingleton::operator[](const std::string& symbol) const {
  auto type = ElementInfo::elementTypeForSymbol(symbol);
  return lookup(type);
}

const ElementDataSingleton& ElementDataSingleton::instance() {
  static ElementDataSingleton instance;
  return instance;
}

const ElementDataSingleton::ElementData& ElementDataSingleton::lookup(const ElementType e) {
  return data().at(e);
}

const std::unordered_map<ElementType, ElementDataSingleton::ElementData>& ElementDataSingleton::data() {
  // clang-format off
  static const std::unordered_map<ElementType, ElementData> map = {
    {ElementType::none,  ElementData("None",    0, 0.0,  0.0, 0.0, -1, 0, 0,  0,  0)},
    {ElementType::H,  ElementData("H",    1, 1.0079,  32, 109, -1, 1, 0,  0,  0)},
    {ElementType::He, ElementData("He",   2, 4.0026,  37, 140, -1, 2, 0,  0,  0)},
    {ElementType::Li, ElementData("Li",   3,  6.941, 130, 182, -1, 1, 0,  0,  0)},
    {ElementType::Be, ElementData("Be",   4, 9.0122,  99, 153, -1, 2, 0,  0,  0)},
    {ElementType::B,  ElementData("B",    5, 10.811,  84, 192, -1, 2, 1,  0,  0)},
    {ElementType::C,  ElementData("C",    6, 12.011,  75, 170, -1, 2, 2,  0,  0)},
    {ElementType::N,  ElementData("N",    7, 14.007,  71, 155, -1, 2, 3,  0,  0)},
    {ElementType::O,  ElementData("O",    8, 15.999,  64, 152, -1, 2, 4,  0,  0)},
    {ElementType::F,  ElementData("F",    9, 18.988,  60, 147, -1, 2, 5,  0,  0)},
    {ElementType::Ne, ElementData("Ne",  10, 20.180,  62, 154, -1, 2, 6,  0,  0)},
    {ElementType::Na, ElementData("Na",  11, 22.990, 160, 227, -1, 1, 0,  0,  0)},
    {ElementType::Mg, ElementData("Mg",  12, 24.305, 140, 173, -1, 2, 0,  0,  0)},
    {ElementType::Al, ElementData("Al",  13, 26.982, 124, 184, -1, 2, 1,  0,  0)},
    {ElementType::Si, ElementData("Si",  14, 28.086, 114, 210, -1, 2, 2,  0,  0)},
    {ElementType::P,  ElementData("P",   15, 30.974, 109, 180, -1, 2, 3,  0,  0)},
    {ElementType::S,  ElementData("S",   16, 32.065, 104, 180, -1, 2, 4,  0,  0)},
    {ElementType::Cl, ElementData("Cl",  17, 35.453, 100, 175, -1, 2, 5,  0,  0)},
    {ElementType::Ar, ElementData("Ar",  18, 39.948, 101, 188, -1, 2, 6,  0,  0)},
    {ElementType::K,  ElementData("K",   19, 39.098, 200, 275, -1, 1, 0,  0,  0)},
    {ElementType::Ca, ElementData("Ca",  20, 40.078, 174, 231, -1, 2, 0,  0,  0)},
    {ElementType::Sc, ElementData("Sc",  21, 44.956, 159, 215, -1, 2, 0,  1,  0)},
    {ElementType::Ti, ElementData("Ti",  22, 47.867, 148, 211, -1, 2, 0,  2,  0)},
    {ElementType::V,  ElementData("V",   23, 50.942, 144, 207, -1, 2, 0,  3,  0)},
    {ElementType::Cr, ElementData("Cr",  24, 51.996, 130, 206, -1, 2, 0,  4,  0)},
    {ElementType::Mn, ElementData("Mn",  25, 54.938, 129, 205, -1, 2, 0,  5,  0)},
    {ElementType::Fe, ElementData("Fe",  26, 55.938, 124, 204, -1, 2, 0,  6,  0)},
    {ElementType::Co, ElementData("Co",  27, 58.933, 118, 200, -1, 2, 0,  7,  0)},
    {ElementType::Ni, ElementData("Ni",  28, 58.693, 117, 197, -1, 2, 0,  8,  0)},
    {ElementType::Cu, ElementData("Cu",  29, 63.546, 122, 196, -1, 2, 0,  9,  0)},
    {ElementType::Zn, ElementData("Zn",  30,  65.38, 120, 201, -1, 2, 0, 10,  0)},
    {ElementType::Ga, ElementData("Ga",  31, 69.723, 123, 187, -1, 2, 1, 10,  0)},
    {ElementType::Ge, ElementData("Ge",  32,  72.64, 120, 211, -1, 2, 2, 10,  0)},
    {ElementType::As, ElementData("As",  33, 74.922, 120, 185, -1, 2, 3, 10,  0)},
    {ElementType::Se, ElementData("Se",  34,  78.96, 118, 190, -1, 2, 4, 10,  0)},
    {ElementType::Br, ElementData("Br",  35, 79.904, 117, 185, -1, 2, 5, 10,  0)},
    {ElementType::Kr, ElementData("Kr",  36, 83.798, 116, 202, -1, 2, 6, 10,  0)},
    {ElementType::Rb, ElementData("Rb",  37, 83.468, 215, 303, -1, 1, 0,  0,  0)},
    {ElementType::Sr, ElementData("Sr",  38,  87.62, 190, 249, -1, 2, 0,  0,  0)},
    {ElementType::Y,  ElementData("Y",   39, 88.906, 176, 232, -1, 2, 0,  1,  0)},
    {ElementType::Zr, ElementData("Zr",  40, 91.224, 164, 223, -1, 2, 0,  2,  0)},
    {ElementType::Nb, ElementData("Nb",  41, 92.906, 156, 218, -1, 2, 0,  3,  0)},
    {ElementType::Mo, ElementData("Mo",  42,  95.96, 146, 217, -1, 2, 0,  4,  0)},
    {ElementType::Tc, ElementData("Tc",  43,  98.91, 138, 216, -1, 2, 0,  5,  0)},
    {ElementType::Ru, ElementData("Ru",  44, 101.07, 136, 213, -1, 2, 0,  6,  0)},
    {ElementType::Rh, ElementData("Rh",  45, 102.91, 134, 210, -1, 2, 0,  7,  0)},
    {ElementType::Pd, ElementData("Pd",  46, 106.42, 130, 210, -1, 2, 0,  8,  0)},
    {ElementType::Ag, ElementData("Ag",  47, 107.87, 136, 211, -1, 2, 0,  9,  0)},
    {ElementType::Cd, ElementData("Cd",  48, 112.41, 140, 218, -1, 2, 0, 10,  0)},
    {ElementType::In, ElementData("In",  49, 114.82, 142, 193, -1, 2, 1, 10,  0)},
    {ElementType::Sn, ElementData("Sn",  50, 118.71, 140, 217, -1, 2, 2, 10,  0)},
    {ElementType::Sb, ElementData("Sb",  51, 121.76, 140, 206, -1, 2, 3, 10,  0)},
    {ElementType::Te, ElementData("Te",  52, 127.60, 137, 206, -1, 2, 4, 10,  0)},
    {ElementType::I,  ElementData("I",   53, 126.90, 136, 198, -1, 2, 5, 10,  0)},
    {ElementType::Xe, ElementData("Xe",  54, 131.29, 136, 216, -1, 2, 6, 10,  0)},
    {ElementType::Cs, ElementData("Cs",  55, 132.91, 238, 343, -1, 1, 0,  0,  0)},
    {ElementType::Ba, ElementData("Ba",  56, 137.33, 206, 268, -1, 2, 0,  0,  0)},
    {ElementType::La, ElementData("La",  57, 138.91, 194, 243, -1, 2, 0,  0,  1)},
    {ElementType::Ce, ElementData("Ce",  58, 140.12, 184, 242, -1, 2, 0,  0,  2)},
    {ElementType::Pr, ElementData("Pr",  59, 140.91, 190, 240, -1, 2, 0,  0,  3)},
    {ElementType::Nd, ElementData("Nd",  60, 144.24, 188, 239, -1, 2, 0,  0,  4)},
    {ElementType::Pm, ElementData("Pm",  61, 146.90, 186, 238, -1, 2, 0,  0,  5)},
    {ElementType::Sm, ElementData("Sm",  62, 150.36, 185, 236, -1, 2, 0,  0,  6)},
    {ElementType::Eu, ElementData("Eu",  63, 151.96, 183, 235, -1, 2, 0,  0,  7)},
    {ElementType::Gd, ElementData("Gd",  64, 157.25, 182, 234, -1, 2, 0,  0,  8)},
    {ElementType::Tb, ElementData("Tb",  65, 158.93, 181, 233, -1, 2, 0,  0,  9)},
    {ElementType::Dy, ElementData("Dy",  66, 162.50, 180, 231, -1, 2, 0,  0, 10)},
    {ElementType::Ho, ElementData("Ho",  67, 164.93, 179, 230, -1, 2, 0,  0, 11)},
    {ElementType::Er, ElementData("Er",  68, 167.26, 177, 229, -1, 2, 0,  0, 12)},
    {ElementType::Tm, ElementData("Tm",  69, 168.93, 177, 227, -1, 2, 0,  0, 13)},
    {ElementType::Yb, ElementData("Yb",  70, 173.05, 178, 226, -1, 2, 0,  0, 14)},
    {ElementType::Lu, ElementData("Lu",  71, 174.97, 174, 224, -1, 2, 0,  1, 14)},
    {ElementType::Hf, ElementData("Hf",  72, 178.49, 164, 223, -1, 2, 0,  2, 14)},
    {ElementType::Ta, ElementData("Ta",  73, 180.95, 158, 222, -1, 2, 0,  3, 14)},
    {ElementType::W,  ElementData("W",   74, 183.84, 150, 218, -1, 2, 0,  4, 14)},
    {ElementType::Re, ElementData("Re",  75, 186.21, 141, 216, -1, 2, 0,  5, 14)},
    {ElementType::Os, ElementData("Os",  76, 190.23, 136, 216, -1, 2, 0,  6, 14)},
    {ElementType::Ir, ElementData("Ir",  77, 192.22, 132, 213, -1, 2, 0,  7, 14)},
    {ElementType::Pt, ElementData("Pt",  78, 195.08, 130, 213, -1, 2, 0,  8, 14)},
    {ElementType::Au, ElementData("Au",  79, 196.97, 130, 214, -1, 2, 0,  9, 14)},
    {ElementType::Hg, ElementData("Hg",  80, 200.59, 132, 223, -1, 2, 0, 10, 14)},
    {ElementType::Tl, ElementData("Tl",  81, 204.38, 144, 196, -1, 2, 1, 10, 14)},
    {ElementType::Pb, ElementData("Pb",  82,  207.2, 145, 202, -1, 2, 2, 10, 14)},
    {ElementType::Bi, ElementData("Bi",  83, 208.98, 150, 207, -1, 2, 3, 10, 14)},
    {ElementType::Po, ElementData("Po",  84, 209.98, 142, 197, -1, 2, 4, 10, 14)},
    {ElementType::At, ElementData("At",  85,    210, 148, 202, -1, 2, 5, 10, 14)},
    {ElementType::Rn, ElementData("Rn",  86,    222, 146, 220, -1, 2, 6, 10, 14)},
    {ElementType::Fr, ElementData("Fr",  87,    223, 242, 348, -1, 1, 0,  0,  0)},
    {ElementType::Ra, ElementData("Ra",  88, 226.03, 211, 283, -1, 2, 0,  0,  0)},
    {ElementType::Ac, ElementData("Ac",  89,    227, 201, 247, -1, 2, 0,  0,  1)},
    {ElementType::Th, ElementData("Th",  90, 232.04, 190, 245, -1, 2, 0,  0,  2)},
    {ElementType::Pa, ElementData("Pa",  91, 231.04, 184, 243, -1, 2, 0,  0,  3)},
    {ElementType::U,  ElementData("U",   92, 238.03, 183, 241, -1, 2, 0,  0,  4)},
    {ElementType::Np, ElementData("Np",  93, 237.05, 180, 239, -1, 2, 0,  0,  5)},
    {ElementType::Pu, ElementData("Pu",  94, 244.10, 180, 243, -1, 2, 0,  0,  6)},
    {ElementType::Am, ElementData("Am",  95, 243.10, 173, 244, -1, 2, 0,  0,  7)},
    {ElementType::Cm, ElementData("Cm",  96, 247.10, 168, 245, -1, 2, 0,  0,  8)},
    {ElementType::Bk, ElementData("Bk",  97, 247.10, 168, 244, -1, 2, 0,  0,  9)},
    {ElementType::Cf, ElementData("Cf",  98, 251.10, 168, 245, -1, 2, 0,  0, 10)},
    {ElementType::Es, ElementData("Es",  99, 254.10, 165, 245, -1, 2, 0,  0, 11)},
    {ElementType::Fm, ElementData("Fm", 100, 257.10, 167, 245, -1, 2, 0,  0, 12)},
    {ElementType::Md, ElementData("Md", 101,    258, 173, 246, -1, 2, 0,  0, 13)},
    {ElementType::No, ElementData("No", 102,    259, 176, 246, -1, 2, 0,  0, 14)},
    {ElementType::Lr, ElementData("Lr", 103,    262, 161, 246, -1, 2, 0,  1, 14)},
    {ElementType::Rf, ElementData("Rf", 104,    261, 157,  -1, -1, 2, 0,  2, 14)},
    {ElementType::Db, ElementData("Db", 105,    262, 149,  -1, -1, 2, 0,  3, 14)},
    {ElementType::Sg, ElementData("Sg", 106,    266, 143,  -1, -1, 2, 0,  4, 14)},
    {ElementType::Bh, ElementData("Bh", 107,    264, 141,  -1, -1, 2, 0,  5, 14)},
    {ElementType::Hs, ElementData("Hs", 108,    277, 134,  -1, -1, 2, 0,  6, 14)},
    {ElementType::Mt, ElementData("Mt", 109,    268, 129,  -1, -1, 2, 0,  7, 14)},
    {ElementType::Ds, ElementData("Ds", 110,    281, 128,  -1, -1, 2, 0,  8, 14)},
    {ElementType::Rg, ElementData("Rg", 111,    280, 121,  -1, -1, 2, 0,  9, 14)},
    {ElementType::Cn, ElementData("Cn", 112,    285, 122,  -1, -1, 2, 0, 10, 14)},
  };
  // clang-format on
  return map;
}

} /* namespace Constants */
} /* namespace Utils */
} /* namespace Scine */
