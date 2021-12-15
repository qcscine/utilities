/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "BSplineTools.h"

namespace Scine {
namespace Utils {

namespace BSplines {

BSplineTools::BSplineTools() = default;

int BSplineTools::findIdxOfLeftOrEqualDomainKnot(const double u, const int p, const Eigen::VectorXd& U) {
  int i = p;

  while ((u >= U(i + 1)) && (i + 1 < U.size() - p)) {
    ++i;
    if (U(i + 1) == U(i)) {
      break; // stop iteration
    }
  }
  return i;
}

int BSplineTools::findIdxOfLeftDomainKnot(const double u, const int p, const Eigen::VectorXd& U) {
  int i = p;
  while ((u > U(i + 1)) && (i + 1 < U.size() - p)) {
    ++i;
  }
  return i;
}

int BSplineTools::findIdxOfRightDomainKnot(const double u, const int p, const Eigen::VectorXd& U) {
  int i = static_cast<int>(U.size()) - 1 - p;
  while ((u < U(i - 1)) && (i - 1 > p)) {
    --i;
  }
  return i;
}

int BSplineTools::findIdxOfRightOrEqualDomainKnot(const double u, const int p, const Eigen::VectorXd& U) {
  int i = static_cast<int>(U.size()) - 1 - p;
  while ((u <= U(i - 1)) && (i - 1 >= p)) {
    --i;
    if (U(i - 1) == U(i)) {
      break; // stop iteration
    }
  }
  return i;
}

double BSplineTools::knotAverage(const int i, const int p, const Eigen::VectorXd& U) {
  double knotSum = 0.0;
  for (int j = i + 1; j <= i + p; ++j) {
    knotSum += U(j);
  }
  return knotSum / double(p);
}

void BSplineTools::normalizeKnotVector(Eigen::VectorXd& knotVector) {
  std::pair<double, double> oldLim = {knotVector.minCoeff(), knotVector.maxCoeff()};

  rescaleKnotVector(knotVector, oldLim, {0, 1});
}

Eigen::VectorXd BSplineTools::normalizedKnotVector(const Eigen::VectorXd& knotVector) {
  std::pair<double, double> oldLim = {knotVector.minCoeff(), knotVector.maxCoeff()};

  return rescaledKnotVector(knotVector, oldLim, {0, 1});
}

double BSplineTools::rescaledKnot(const double knot, const std::pair<double, double> oldLim,
                                  const std::pair<double, double> newLim) {
  double knotCopy = knot;
  knotCopy -= oldLim.first;
  knotCopy /= (oldLim.second - oldLim.first);
  knotCopy *= (newLim.second - newLim.first);
  knotCopy += newLim.first;
  return knotCopy;
}

Eigen::VectorXd BSplineTools::rescaledKnotVector(const Eigen::VectorXd& knotVector, const std::pair<double, double> oldLim,
                                                 const std::pair<double, double> newLim) {
  auto Ucopy = knotVector;
  rescaleKnotVector(Ucopy, oldLim, newLim);
  return Ucopy;
}

void BSplineTools::rescaleKnotVector(Eigen::VectorXd& knotVector, const std::pair<double, double> oldLim,
                                     const std::pair<double, double> newLim) {
  knotVector.array() -= oldLim.first;
  knotVector /= (oldLim.second - oldLim.first);
  knotVector *= (newLim.second - newLim.first);
  knotVector.array() += newLim.first;
}

// finite difference matrix penalizing BSplines
int BSplineTools::differenceOperator(int i, int j, int kappa) {
  if (kappa > 1) {
    return differenceOperator(i + 1, j, kappa - 1) - differenceOperator(i, j, kappa - 1);
  }
  if (kappa == 1) {
    if (i + 1 == j) {
      return 1;
    }
    if (i == j) {
      return -1;
    }

    return 0;
  }
  return 0;
}

} // namespace BSplines
} // namespace Utils
} // namespace Scine
