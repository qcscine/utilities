/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef BSPLINES_GENERATOR_H
#define BSPLINES_GENERATOR_H

#include "BSpline.h"
#include <Eigen/Core>

namespace Scine {
namespace Utils {

namespace BSplines {

/*!
 * Basis class for BSpline generators.
 * NB: We do not take a "const Eigen::MatrixXd&" since this would create a temporary when using it with Eigen::VectorXd.
 */
class Generator {
 public:
  Generator(const Eigen::Ref<const Eigen::MatrixXd>& dataPoints, int splineDegree);
  Generator(const Eigen::Ref<const Eigen::MatrixXd>& dataPoints, int numberOfControlPoints, int splineDegree);
  virtual ~Generator() = default;

  BSpline generateBSpline();

 protected:
  const Eigen::Ref<const Eigen::MatrixXd>& dataPoints_;
  int p_, dim_, m_, n_;
  Eigen::VectorXd U_, uBar_;
  Eigen::MatrixXd controlPoints_;

 private:
  virtual Eigen::VectorXd generateKnotVector() = 0;
  virtual Eigen::MatrixXd generateControlPoints() = 0;
};

} // namespace BSplines

} // namespace Utils
} // namespace Scine
#endif // BSPLINES_GENERATOR_H
