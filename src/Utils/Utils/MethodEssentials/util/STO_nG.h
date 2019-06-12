/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_STO_NG_H
#define UTILS_STO_NG_H

#include <Utils/MethodEssentials/util/GTOExpansion.h>
#include <array>
#include <vector>

namespace Scine {
namespace Utils {

/*!
 * The class STO_nG delivers the coefficients of a STO-nG expansion for a
 * given STO.
 * Expansion coefficients reference:
 * 1s-5g: Robert F. Stewart, Small Gaussian Expansions of Slater‐Type Orbitals, J. Chem. Phys. 52, 431 (1970).
 * >= 6s: Prof. Dr. Jürg Hutter.
 */

class STO_nG {
 public:
  using parameterPair = std::pair<double, double>;

  static std::vector<parameterPair> get(unsigned int N, unsigned int n, unsigned int l, double exponent = 1.0);
  static GTOExpansion getGTOExpansion(unsigned int N, unsigned int n, unsigned int l, double exponent = 1.0);

 private:
  static void getValues(std::array<double, 6>& a, std::array<double, 6>& c, int N, int n, int l);
  static void sto1s(std::array<double, 6>& a, std::array<double, 6>& c, int N);
  static void sto2s(std::array<double, 6>& a, std::array<double, 6>& c, int N);
  static void sto2p(std::array<double, 6>& a, std::array<double, 6>& c, int N);
  static void sto3s(std::array<double, 6>& a, std::array<double, 6>& c, int N);
  static void sto3p(std::array<double, 6>& a, std::array<double, 6>& c, int N);
  static void sto3d(std::array<double, 6>& a, std::array<double, 6>& c, int N);
  static void sto4s(std::array<double, 6>& a, std::array<double, 6>& c, int N);
  static void sto4p(std::array<double, 6>& a, std::array<double, 6>& c, int N);
  static void sto4d(std::array<double, 6>& a, std::array<double, 6>& c, int N);
  static void sto4f(std::array<double, 6>& a, std::array<double, 6>& c, int N);
  static void sto5s(std::array<double, 6>& a, std::array<double, 6>& c, int N);
  static void sto5p(std::array<double, 6>& a, std::array<double, 6>& c, int N);
  static void sto5d(std::array<double, 6>& a, std::array<double, 6>& c, int N);
  static void sto5f(std::array<double, 6>& a, std::array<double, 6>& c, int N);
  static void sto5g(std::array<double, 6>& a, std::array<double, 6>& c, int N);
  static void sto6s(std::array<double, 6>& a, std::array<double, 6>& c, int N);
  static void sto6p(std::array<double, 6>& a, std::array<double, 6>& c, int N);
  static void sto6d(std::array<double, 6>& a, std::array<double, 6>& c, int N);
  static void sto6f(std::array<double, 6>& a, std::array<double, 6>& c, int N);
  static void sto6g(std::array<double, 6>& a, std::array<double, 6>& c, int N);
  static void sto6h(std::array<double, 6>& a, std::array<double, 6>& c, int N);
  static void sto7s(std::array<double, 6>& a, std::array<double, 6>& c, int N);
  static void sto7p(std::array<double, 6>& a, std::array<double, 6>& c, int N);
  static void sto7d(std::array<double, 6>& a, std::array<double, 6>& c, int N);
  static void sto7f(std::array<double, 6>& a, std::array<double, 6>& c, int N);
  static void sto7g(std::array<double, 6>& a, std::array<double, 6>& c, int N);
  static void sto7h(std::array<double, 6>& a, std::array<double, 6>& c, int N);
  static void sto7i(std::array<double, 6>& a, std::array<double, 6>& c, int N);
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_STO_NG_H
