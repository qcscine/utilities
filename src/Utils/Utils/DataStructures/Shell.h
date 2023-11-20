/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILSOS_SHELL_H
#define UTILSOS_SHELL_H

#include <Utils/Typenames.h>

namespace Scine {
namespace Utils {
namespace Integrals {

/** Generally contracted Solid-Harmonic/Cartesian Gaussian Shell class implementation
 *  Using the standard libint shell ordering.
 *
 *  For real-solid harmomnics:
 *  m=-l,-l+1,...,0,...,l-1,l
 *  See Helgaker Molecular Electronic Structure Theory. Eqn. 6.4.19-21
 *
 *  For cartesian Gaussians:
 *  (See Libint Programmers manual)
 *  for (int i = 0; i <= l; ++i) {
 *    int lx = l - i; // exponent of x
 *    for (int j = 0; j <= i; ++j) {
 *      int ly = i - j; // exponent of y
 *      int lz = j;     // exponent of z
 *      // 2.
 *      // phi(r) = S * x^lx y^ly z^lz
 *      dens(d) = S * cartesian(lx, ly, lz, dx, dy, dz);
 *      ++d;
 *    }
 *  }
 *
 *  Note: The difference to Libint is that we provide only one max angular momentum per shell.
 *
 */
class Shell {
 private:
  std::vector<double> _vecAlpha;    // vector of Gaussian widths
  std::vector<double> _vecCoeffs;   // vector of contraction coefficients
  std::vector<double> _maxLnCoeffs; // max log of contraction coefficients
  Displacement _shift;              // Gaussian shift -- equal for all GTFs in Shell
  std::size_t _l;                   // angular momentum
  bool _pureSolid;                  // true if solid harmonics are used
  std::size_t _contrLegnth;         // Number of contracted Gaussian type functions.

 public:
  Shell();
  ~Shell();

  /**
   * @brief Constructor of a shell.
   * The Shell is defined as a contracted Gaussian type orbital. Either with solid harmonics or cartesians prefactors.
   * Note that, opposed to Libint, our Shell is given only by a single angular momentum.
   * Therefore, e.g., a single Libint STO-3G s+p shell object, is represented by two shell objects here.
   *
   * @note Do not normalize the contraction coefficients and let Libint handle this!
   *       Special care must be taken, when Shell objects are converted to Libint, since Libint automatically
   *       renormalizes the coefficients, which will lead to wrong results.
   *
   *
   * @param v_alpha
   * @param v_coeffs
   * @param shift
   * @param max_l
   * @param pure_solid
   */
  Shell(std::vector<double> v_alpha, std::vector<double> v_coeffs, Displacement shift, std::size_t max_l, bool pure_solid);

  /**
   * @brief Returns the size of the Shell.
   * It takes into account whether it is a solid harmonic or cartesian shell.
   */
  inline auto size() const -> std::size_t {
    return _pureSolid ? (2 * _l + 1) : _cartesianSize();
  }
  /**
   * @brief Returns the number of primitives in shell.
   */
  inline auto nprim() const -> std::size_t {
    return _vecAlpha.size();
  }
  /**
   * @brief Moves the shell to ``new_shit``
   * @param new_shift
   */
  inline auto move(const Displacement& new_shift) -> void {
    this->_shift = new_shift;
  }
  /**
   * @brief Overwrite the alpha parameter vector.
   * @note Must be the same length as the contraction.
   * @param new_alpha
   */
  inline auto setAlpha(const std::vector<double>& new_alpha) -> void {
    this->_vecAlpha = new_alpha;
  }
  /**
   * @brief Overwrite the coefficient vector.
   * @note Must be the same length as the contraction.
   * @param new_coeffs
   */
  inline auto setCoefficients(const std::vector<double>& new_coeffs) -> void {
    this->_vecCoeffs = new_coeffs;
  }

  inline const std::vector<double>& getVecAlpha() const {
    return _vecAlpha;
  }

  inline const std::vector<double>& getVecCoeffs() const {
    return _vecCoeffs;
  }

  inline const std::vector<double>& getMaxLnCoeffs() const {
    return _maxLnCoeffs;
  }

  inline const Displacement& getShift() const {
    return _shift;
  }

  inline size_t l() const {
    return _l;
  }

  inline bool isPureSolid() const {
    return _pureSolid;
  }

  inline size_t getContrLegnth() const {
    return _contrLegnth;
  }

 private:
  inline auto _cartesianSize() const -> std::size_t {
    return (_l + 1) * (_l + 2) / 2;
  }

  /**
   * @brief Evaluate max log coefficient vector
   */
  inline void _calcMaxLogCoeffs() {
    _maxLnCoeffs.resize(_contrLegnth);
    for (auto p = 0UL; p < _contrLegnth; ++p) {
      double max_ln_c = -std::numeric_limits<double>::max();
      max_ln_c = std::max(max_ln_c, std::log(std::abs(_vecCoeffs[p])));
      _maxLnCoeffs[p] = max_ln_c;
    }
  }
};

} // namespace Integrals
} // namespace Utils
} // namespace Scine

#endif // UTILSOS_SHELL_H
