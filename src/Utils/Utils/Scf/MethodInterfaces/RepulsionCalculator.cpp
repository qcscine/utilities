/**
 * @file RepulsionCalculator.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "RepulsionCalculator.h"
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/Math/AutomaticDifferentiation/MethodsHelpers.h>
#if defined(_OPENMP)
#  include <omp.h>
#endif

namespace Scine {
namespace Utils {
using namespace AutomaticDifferentiation;

RepulsionCalculator::RepulsionCalculator(const Utils::ElementTypeCollection& elements, const Utils::PositionCollection& positions)
  : elements_(elements), positions_(positions), repulsionEnergy_(0.0) {
}

void RepulsionCalculator::initialize() {
  // Initialize one to make sure the singleton is generated.
  double repulsionConstant = ElementInfo::Z(elements_[0]);
}

void RepulsionCalculator::calculateRepulsion(Utils::DerivativeOrder order) {
  repulsionEnergy_ = 0;
  repulsionGradients_ = GradientCollection::Zero(elements_.size(), 3);

  // prepare charges
  std::vector<double> charges(elements_.size());
  for (int i = 0; i < elements_.size(); ++i)
    charges[i] = ElementInfo::Z(elements_[i]);

  // prepare hessian
  if (order == DerivativeOrder::Two) {
    AutomaticDifferentiation::Second3D dummy;
    for (int i = 0; i < elements_.size(); ++i) {
      for (int j = i + 1; j < elements_.size(); ++j) {
        repulsionHessian_.insert({{i, j}, dummy});
      }
    }
  }

#if defined(_OPENMP)
  const int nThreads = omp_get_max_threads();
#else
  const int nThreads = 1;
#endif

  std::vector<double> repulsionBuffer(nThreads, 0.0);
  std::vector<GradientCollection> gradientBuffer(nThreads, repulsionGradients_);

  // parallel evaluation
  for (int i = 0; i < elements_.size(); ++i) {
#pragma omp parallel for schedule(static)
    for (int j = i + 1; j < elements_.size(); ++j) {
#if defined(_OPENMP)
      const int thread = omp_get_thread_num();
#else
      const int thread = 0;
#endif
      Eigen::RowVector3d distanceVector = positions_.row(j) - positions_.row(i);
      const double distance = distanceVector.norm();
      const double repulsionConstant = charges[i] * charges[j];

      if (order == DerivativeOrder::Zero) {
        repulsionBuffer[thread] += calculatePairwiseCoreRepulsion<DerivativeOrder::Zero>(distance, repulsionConstant);
      }
      else if (order == DerivativeOrder::One) {
        First1D gradientData = calculatePairwiseCoreRepulsion<Utils::DerivativeOrder::One>(distance, repulsionConstant);
        repulsionBuffer[thread] += gradientData.value();
        Gradient pairwiseGradient =
            Gradient(get3Dfrom1D<DerivativeOrder::One>(gradientData, distanceVector).derivatives().transpose());
        addDerivativeToContainer<Derivative::First>(gradientBuffer[thread], i, j, pairwiseGradient);
      }
      else if (order == DerivativeOrder::Two) {
        auto pairwiseHessian = get3Dfrom1D<DerivativeOrder::Two>(
            calculatePairwiseCoreRepulsion<DerivativeOrder::Two>(distance, repulsionConstant), distanceVector);
        repulsionBuffer[thread] += pairwiseHessian.value();
        addDerivativeToContainer<Derivative::First>(gradientBuffer[thread], i, j, Gradient(pairwiseHessian.deriv()));
        repulsionHessian_.at({i, j}) = pairwiseHessian;
      }
    }
  }
  // reduction of parallel region buffers
  for (int i = 0; i < nThreads; ++i) {
    repulsionEnergy_ += repulsionBuffer[i];
    repulsionGradients_ += gradientBuffer[i];
  }
}

void RepulsionCalculator::addRepulsionDerivatives(DerivativeContainerType<Utils::Derivative::First>& derivatives) const {
  derivatives += repulsionGradients_;
}

void RepulsionCalculator::addRepulsionDerivatives(DerivativeContainerType<Utils::Derivative::SecondAtomic>& derivatives) const {
  for (int i = 0; i < elements_.size(); ++i) {
    for (int j = i + 1; j < elements_.size(); ++j) {
      addDerivativeToContainer<Derivative::SecondAtomic>(derivatives, i, j, repulsionHessian_.at({i, j}));
    }
  }
}

void RepulsionCalculator::addRepulsionDerivatives(DerivativeContainerType<Utils::Derivative::SecondFull>& derivatives) const {
  for (int i = 0; i < elements_.size(); ++i) {
    for (int j = i + 1; j < elements_.size(); ++j) {
      addDerivativeToContainer<Derivative::SecondFull>(derivatives, i, j, repulsionHessian_.at({i, j}));
    }
  }
}

double RepulsionCalculator::getRepulsionEnergy() const {
  return repulsionEnergy_;
}

} // namespace Utils
} // namespace Scine
