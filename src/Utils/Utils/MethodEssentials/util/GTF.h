/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_GTF_H
#define UTILS_GTF_H

namespace Scine {
namespace Utils {

/*!
 * Class that contains one of the gaussian functions
 * of a ,f.i., STO-nG expansion
 */

class GTF {
 public:
  GTF() = default;

  GTF(int l, double a, double c) {
    set(l, a, c);
  }
  void set(int l, double a, double c);
  double getExponent() const {
    return exponent_;
  }
  double getCoefficient() const {
    return coefficient_;
  }
  double getNormalizedCoefficient() const {
    return normalizedCoefficient_;
  }

 private:
  void normalizeCoefficient(int l, double c);

  double exponent_;
  double coefficient_;
  double normalizedCoefficient_;
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_GTF_H
