/**
 * @file ElementTypes.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_ELEMENTTYPES_H_
#define UTILS_ELEMENTTYPES_H_

namespace Scine {
namespace Utils {

// clang-format off
/**
 * @brief Enum class defining all elements
 */
enum class ElementType {
  none = 0,
  H,                                                                                                                          He,
  Li, Be,                                                                                                 B,  C,  N,  O,  F,  Ne,
  Na, Mg,                                                                                                 Al, Si, P,  S,  Cl, Ar,
  K,  Ca, Sc,                                                         Ti, V,  Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr,
  Rb, Sr, Y,                                                          Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I,  Xe,
  Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta, W,  Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn,
  Fr, Ra, Ac, Th, Pa, U,  Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn
};
// clang-format on
} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_ELEMENTTYPES_H_
