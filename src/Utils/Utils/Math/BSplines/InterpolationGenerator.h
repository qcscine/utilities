/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef BSPLINES_INTERPOLATIONGENERATOR_H
#define BSPLINES_INTERPOLATIONGENERATOR_H

#include "Generator.h"
#include <Eigen/SVD>

namespace Scine {
namespace Utils {

namespace BSplines {

/*!
 * Generates a B-spline curve interpolating the data points.
 * (see the NURBS book by Piegl 1997)
 */
class InterpolationGenerator : public Generator {
 public:
  /*! @param dataPoints number rows = number of points, number of cols = dimension */
  explicit InterpolationGenerator(const Eigen::Ref<const Eigen::MatrixXd>& dataPoints, unsigned degree = 3,
                                  bool uniformKnotVector = false);

 private:
  Eigen::VectorXd generateKnotVector() override;
  Eigen::MatrixXd generateControlPoints() override;
  Eigen::MatrixXd generateControlPointMatrix();
  void calculateCoefficientMatrix();
  void initializeSolver();

  bool uniformKnotVector_;
  Eigen::MatrixXd Nmat_;
  Eigen::JacobiSVD<Eigen::MatrixXd> svdOfNmat_;
};

} // namespace BSplines

} // namespace Utils
} // namespace Scine
#endif // BSPLINES_INTERPOLATIONGENERATOR_H
