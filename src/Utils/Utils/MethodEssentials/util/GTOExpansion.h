/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_GTOEXPANSION_H
#define UTILS_GTOEXPANSION_H

#include <Utils/MethodEssentials/util/GTF.h>
#include <vector>

namespace Scine {
namespace Utils {

/*!
 * The class GTOExpansion is the container for the coefficients of a
 * STO-nG expansion.
 */

class GTOExpansion {
 public:
  explicit GTOExpansion(int N = 0) : par_(static_cast<unsigned>(N)) {
    setAngularMomentum(0);
  }

  void setSize(int N) {
    par_.resize(N);
  }
  int size() const {
    return static_cast<int>(par_.size());
  }

  void setAngularMomentum(int l) {
    angularMomentum_ = l;
    nAOs_ = (l == 0) ? 1 : (l == 1) ? 3 : 5;
  }
  int angularMomentum() const {
    return angularMomentum_;
  }
  int nAOs() const {
    return nAOs_;
  }

  void setGTF(int index, const GTF& p) {
    par_[index] = p;
  }
  void setGTF(int index, double a, double c) {
    par_[index] = GTF(angularMomentum_, a, c);
  }
  const GTF& getGTF(int index) const {
    return par_[index];
  }

  // void setExponent(int index, double v) { par_[index].first = v; }
  double getExponent(int index) const {
    return par_[index].getExponent();
  }
  // void setCoefficient(int index, double v) { par_[index].second = v; }
  double getCoefficient(int index) const {
    return par_[index].getCoefficient();
  }
  double getNormalizedCoefficient(int index) const {
    return par_[index].getNormalizedCoefficient();
  }

 private:
  std::vector<GTF> par_;
  int angularMomentum_ = 0;
  int nAOs_ = 0; // Number of AOs sharing those parameters: 1 for s, 3 for p, 5 for d
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_GTOEXPANSION_H
