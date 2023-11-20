/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
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

struct Gtf {
  // Empty constructor
  Gtf() = default;
  // Constructor setting exponent, coefficient and normalized coefficient
  Gtf(int l, double expo, double coef);

  double exponent = 0;
  double coefficient = 0;
  double normalizedCoefficient = 0;

  // Sets the normalized coefficient from exponent and coefficient
  void setNormalized(int l);
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_GTF_H
