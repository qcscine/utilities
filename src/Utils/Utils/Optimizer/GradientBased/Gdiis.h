/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_GDIIS_H_
#define UTILS_GDIIS_H_

#include <Eigen/Dense>
#include <array>

namespace Scine {
namespace Utils {

/**
 * @brief An implementation of the GDIIS optimization acceleration algorithm.
 */
class Gdiis {
 public:
  /**
   * @brief Construct a new GDIIS.
   * @param hessianInverse A reference to the Hessian inverse to be used for the
   *                       generation of an actual step. The content of the
   *                       reference can continuosly be updated using e.g. the
   *                       BFGS scheme or it can also just remain the identity.
   * @param maxm The maximum number of old steps to be stored and to be
   *             extrapolated from. Upon reaching the maximum number of points,
   *             the oldest one will be replaced with the newly given one.
   */
  Gdiis(const Eigen::MatrixXd& hessianInverse, unsigned int maxm)
    : _invH(hessianInverse),
      _maxm(maxm),
      _nParam(hessianInverse.cols()),
      _parmeters(hessianInverse.cols(), maxm),
      _gradients(hessianInverse.cols(), maxm),
      _errors(hessianInverse.cols(), maxm){};
  /**
   * @brief Store data into the GDIIS without extrapolation.
   * @param parameters The current parameters.
   * @param gradients  The current gradients (for the current parameters).
   */
  void store(Eigen::VectorXd& parameters, Eigen::VectorXd& gradients) {
    unsigned int current = _cycle % _maxm;
    _parmeters.col(current) = parameters;
    _gradients.col(current) = gradients;
    _errors.col(current) = -_invH * gradients;
    _cycle++;
  }
  /**
   * @brief Resets the storage of the GDIIS.
   */
  void flush() {
    _cycle = 0;
  }
  /**
   * @brief Stores the new data and extrapolates to optimum parameters using all
   *        stored data.
   * @param parameters The current parameters.
   * @param gradients  The current gradients (for the current parameters).
   * @return Eigen::VectorXd The extrapolated best parameters.
   */
  Eigen::VectorXd update(Eigen::VectorXd& parameters, Eigen::VectorXd& gradients) {
    unsigned int current = _cycle % _maxm;
    Eigen::VectorXd ref = -_invH * gradients;
    _parmeters.col(current) = parameters;
    _gradients.col(current) = gradients;
    _errors.col(current) = ref;

    unsigned int n = _maxm;
    if (_cycle < _maxm - 1)
      n = _cycle;

    Eigen::VectorXd extrapolation;
    if (_cycle == 0) {
      _cycle++;
      return parameters + ref;
    }
    else {
      _cycle++;
      Eigen::MatrixXd B(n, n);
      for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = 0; j < i + 1; j++) {
          const double tmp = _errors.col(i).dot(_errors.col(j));
          B(i, j) = tmp;
          B(j, i) = tmp;
        }
      }
      // Rescale B such that the smallest error vector has atleast a norm of 1
      double minval = B.diagonal().minCoeff();
      if (minval < 1.0) {
        B *= (1.0 / minval);
      }

      // Generate Coefficients
      Eigen::VectorXd rhs(n);
      rhs.fill(1.0);
      Eigen::VectorXd coefficients(n);
      coefficients = B.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(rhs);
      // Check if generated coefficients are stable
      if (coefficients.array().abs().maxCoeff() > 1.0e+8) {
        return parameters + ref;
      }

      // Extrapolate
      coefficients /= coefficients.sum();
      Eigen::VectorXd tmpGrad = Eigen::VectorXd::Zero(_nParam);
      Eigen::VectorXd tmpParam = Eigen::VectorXd::Zero(_nParam);
      for (unsigned int i = 0; i < n; i++) {
        tmpGrad.noalias() += coefficients[i] * _gradients.col(i);
        tmpParam.noalias() += coefficients[i] * _parmeters.col(i);
      }
      extrapolation = tmpParam - _invH * tmpGrad;

      // Check if the extrapolated parameters are actually sensible
      auto step = extrapolation - _parmeters.col(current);
      if (step.norm() > 5.0 * _errors.col(current).norm()) {
        extrapolation += (5.0 / step.norm() - 1.0) * step;
      }
      if (coefficients.array().abs().sum() > 30.0) {
        return parameters + ref;
      }
      const double cosAngle = step.dot(ref) / (step.norm() * ref.norm());
      constexpr std::array<double, 10> cutoffs{{1.00, 1.00, 0.97, 0.84, 0.71, 0.67, 0.62, 0.56, 0.49, 0.41}};
      const double cutoff = (n < 10) ? cutoffs[n] : 0.0;
      if (cosAngle < 0.0) {
        _cycle = 0;
        return parameters + ref;
      }
      else if (cosAngle < cutoff) {
        return parameters + ref;
      }
      return extrapolation;
    }
  }

 private:
  const Eigen::MatrixXd& _invH;
  const unsigned int _maxm;
  const unsigned int _nParam;
  unsigned int _cycle = 0;
  Eigen::MatrixXd _parmeters;
  Eigen::MatrixXd _gradients;
  Eigen::MatrixXd _errors;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_GDIIS_H_
