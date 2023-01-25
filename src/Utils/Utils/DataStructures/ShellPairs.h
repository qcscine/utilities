/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILSOS_SHELLPAIRS_H
#define UTILSOS_SHELLPAIRS_H

#include "Shell.h"
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <vector>

namespace Scine {
namespace Utils {
namespace Integrals {

/*
 * Note: The ShellPairs must be evaluated with an integral evaluation routine of choice.
 * They may be generated via the libint interface using `generateShellPairs`.
 */

/* ShellPairType pre-computes shell-pair data, primitive pairs are screened to finite precision
 * Copied from Libint
 */
struct ShellPairType {
  struct PrimPairData {
    double P[3]; //!< \f$ (\alpha_1 \vec{A} + \alpha_2 \vec{B})/(\alpha_1 + \alpha_2) \f$
    double K;
    double one_over_gamma;
    double scr;
    int p1;
    int p2;
  };

  std::vector<PrimPairData> primpairs;
  double AB[3];

  ShellPairType() : primpairs() {
    for (int i = 0; i != 3; ++i)
      AB[i] = 0.;
  }

  ShellPairType(size_t max_nprim) : primpairs() {
    primpairs.reserve(max_nprim * max_nprim);
    for (int i = 0; i != 3; ++i)
      AB[i] = 0.;
  }
  template<typename Real>
  ShellPairType(const Shell& s1, const Shell& s2, Real ln_prec) {
    init(s1, s2, ln_prec);
  }

  void resize(std::size_t max_nprim) {
    const auto max_nprim2 = max_nprim * max_nprim;
    if (max_nprim * max_nprim > primpairs.size())
      primpairs.resize(max_nprim2);
  }

  /// initializes "expensive" primitive pair data; a pair of primitives with exponents \f$ \{\alpha_a,\alpha_b\} \f$
  /// located at \f$ \{ \vec{A},\vec{B} \} \f$ whose max coefficients in contractions are \f$ \{ \max{|c_a|} ,
  /// \max{|c_b|}
  /// \} \f$ is screened-out (omitted) if \f$ \exp(-|\vec{A}-\vec{B}|^2 \alpha_a * \alpha_b / (\alpha_a + \alpha_b))
  /// \max{|c_a|} \max{|c_b|} \leq \epsilon \f$ where \f$ \epsilon \f$ is the desired precision of the integrals.
  template<typename Real>
  void init(const Shell& s1, const Shell& s2, Real ln_prec) {
    primpairs.clear();

    const auto& A = s1.getShift();
    const auto& B = s2.getShift();
    double AB2 = 0.;
    for (int i = 0; i != 3; ++i) {
      AB[i] = A[i] - B[i];
      AB2 += AB[i] * AB[i];
    }

    const auto& veca1 = s1.getVecAlpha();
    const auto& veca2 = s2.getVecAlpha();

    size_t c = 0;
    for (size_t p1 = 0; p1 != s1.getContrLegnth(); ++p1) {
      for (size_t p2 = 0; p2 != s2.getContrLegnth(); ++p2) {
        const auto& a1 = veca1[p1];
        const auto& a2 = veca2[p2];
        const auto gamma = a1 + a2;
        const auto oogamma = 1 / gamma;

        const auto rho = a1 * a2 * oogamma;
        const auto minus_rho_times_AB2 = -rho * AB2;
        const auto screen_fac = minus_rho_times_AB2 + s1.getMaxLnCoeffs()[p1] + s2.getMaxLnCoeffs()[p2];
        if (screen_fac < ln_prec) {
          continue;
        }
        primpairs.resize(c + 1);
        PrimPairData& p = primpairs[c];
        p.scr = screen_fac;
        p.p1 = p1;
        p.p2 = p2;
        constexpr decltype(rho) sqrt_two_times_M_PI_to_1pt25 = 5.9149671727956128778; // \sqrt{2} (\pi)^{5/4}
        p.K = sqrt_two_times_M_PI_to_1pt25 * exp(minus_rho_times_AB2) * oogamma;
        // p.K = exp(minus_rho_times_AB2) * oogamma;
        if (AB2 == 0.) { // this buys a bit more precision
          p.P[0] = A[0];
          p.P[1] = A[1];
          p.P[2] = A[2];
        }
        else {
          p.P[0] = (a1 * A[0] + a2 * B[0]) * oogamma;
          p.P[1] = (a1 * A[1] + a2 * B[1]) * oogamma;
          p.P[2] = (a1 * A[2] + a2 * B[2]) * oogamma;
        }
        p.one_over_gamma = oogamma;

        ++c;
      }
    }
  }
};

struct ShellPairData {
  size_t secondShellIndex;
  std::shared_ptr<ShellPairType> precomputedShellPair;
  double cauchySchwarzFactor = 0.0;
};

/**
 * @class ShellPair @file ShellPairs.h
 * @brief Class containing the list of all Shells relevant for a shell.
 * This class contains 3 vectors:
 * interactions:          contains the indices of all the shells that form a non negligible
 *                        shell pair with this shell.
 * precomputedShellPairs: contains the actual contracted shell pair as constructed by libint2.
 *                        This is needed by the libint2::compute2 function.
 * cauchySchwarzFactors:  contains the Cauchy-Schwarz factors of a given shell pair, \f$ \sqrt{(\mu\nu|\mu\nu)} \f$ .
 */
class ShellPair {
 public:
  using iterator = std::vector<ShellPairData>::iterator;
  using const_iterator = std::vector<ShellPairData>::const_iterator;

  ShellPair() = default;
  /**
   * @brief Constructor. Assigns a const reference to the shell this list belongs to.
   * @param shell Const ref to a shell to which a list of interacting shell has to be constructed.
   * @param calculateCauchySchwarzFactor determines whether the maximal \f$\sqrt{\left(ij|ij\right)}\f$ coulomb
   *                                     integral is calculated for the \f$i\f$, \f$j\f$ shell-pair.
   */
  explicit ShellPair(const Shell& shell, bool calculateCauchySchwarzFactor = true);
  ~ShellPair() = default;
  ShellPair(const ShellPair& rhs) = delete;
  ShellPair(ShellPair&& rhs) = default;
  ShellPair& operator=(const ShellPair& rhs) = delete;
  ShellPair& operator=(ShellPair&& rhs) = default;

  const Shell& getShell() const;

  iterator begin();
  const_iterator begin() const;
  iterator end();
  const_iterator end() const;
  int size() const;
  ShellPairData& front();
  const ShellPairData& front() const;
  ShellPairData& operator[](int index);
  const ShellPairData& operator[](int index) const;
  ShellPairData& at(int index);
  const ShellPairData& at(int index) const;
  template<typename T>
  void emplace_back(T&& newElement) {
    shellPairDataVector_.emplace_back(std::forward<T>(newElement));
  }

 private:
  std::vector<ShellPairData> shellPairDataVector_;
  bool calculateCauchySchwarzFactor_{true};
  std::unique_ptr<Shell> shell_;
};

/**
 * @class ShellPairs @file ShellPairs.h
 * @brief Typedef representing a list of ShellPair.
 */
class ShellPairs {
 public:
  using iterator = std::vector<ShellPair>::iterator;
  using const_iterator = std::vector<ShellPair>::const_iterator;
  ShellPairs() = default;
  ~ShellPairs();
  bool hasCauchySchwarzFactor() const {
    return hasCauchySchwarzFactors_;
  }

  void setCauchySchwarzFactor(bool hasCauchySchwarz) {
    hasCauchySchwarzFactors_ = hasCauchySchwarz;
  }

  iterator begin();
  const_iterator begin() const;
  iterator end();
  const_iterator end() const;
  int size() const;
  bool empty() const;
  void resize(int newSize);
  ShellPair& operator[](int index);
  const ShellPair& operator[](int index) const;
  ShellPair& at(int index);
  const ShellPair& at(int index) const;
  template<typename T>
  void emplace_back(T&& newElement) {
    shellPairs_.emplace_back(std::forward<T>(newElement));
  }

 private:
  bool hasCauchySchwarzFactors_{false};
  std::vector<ShellPair> shellPairs_;
};

} // namespace Integrals
} // namespace Utils
} // namespace Scine

#endif // UTILSOS_SHELLPAIRS_H
