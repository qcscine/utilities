/**
 * @file NumericalHessianCalculator.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "NumericalHessianCalculator.h"
#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/CalculatorBasics/StatesHandler.h>

namespace Scine {
namespace Utils {

#pragma omp declare reduction(+ : Eigen::MatrixXd : omp_out += omp_in) initializer(omp_priv = Eigen::MatrixXd::Zero(omp_orig.rows(), omp_orig.cols()))
#pragma omp declare reduction(+ : DipoleGradient : omp_out += omp_in) initializer(omp_priv = DipoleGradient::Zero(omp_orig.rows(), omp_orig.cols()))

NumericalHessianCalculator::NumericalHessianCalculator(Core::Calculator& calculator) : calculator_(calculator) {
}

void NumericalHessianCalculator::requiredDipoleGradient(bool dipoleGradient) {
  dipoleGradient_ = dipoleGradient;
}

Results NumericalHessianCalculator::calculate(double delta) {
  PropertyList originalProperties = calculator_.getRequiredProperties();
  auto result = calculateFromGradientDifferences(delta);
  calculator_.setRequiredProperties(originalProperties);
  return result;
}

Results NumericalHessianCalculator::calculate(const std::vector<int>& indices, double delta) {
  PropertyList originalProperties = calculator_.getRequiredProperties();
  auto result = calculateFromGradientDifferences(delta, indices);
  calculator_.setRequiredProperties(originalProperties);
  return result;
}

HessianMatrix NumericalHessianCalculator::calculateFromEnergyDifferences(double delta) {
  calculator_.setRequiredProperties(Utils::Property::Energy);
  auto positions = calculator_.getPositions();

  auto nDimensions = 3 * positions.rows();
  Eigen::MatrixXd m = Eigen::MatrixXd::Zero(nDimensions, nDimensions);

  for (int i = 0; i < nDimensions; ++i) {
    m(i, i) = hessianElementSameFromEnergy(i, positions, delta);
    for (int j = 0; j < i; ++j) {
      m(i, j) = hessianElementDifferentFromEnergy(i, j, positions, delta);
      m(j, i) = m(i, j);
    }
  }

  calculator_.modifyPositions(positions);
  return HessianMatrix{m};
}

double NumericalHessianCalculator::hessianElementSameFromEnergy(int i, const PositionCollection& referencePositions,
                                                                double delta) {
  auto modifiedPositions = referencePositions;

  auto atomI = i / 3;
  auto compI = i % 3;

  calculator_.modifyPositions(modifiedPositions);
  auto r = calculator_.calculate("");
  if (!r.get<Property::SuccessfulCalculation>())
    throw std::runtime_error("Energy calculation in numerical Hessian evaluation failed.");
  auto E = r.get<Property::Energy>();

  modifiedPositions.row(atomI)(compI) = referencePositions.row(atomI)(compI) - delta;
  calculator_.modifyPositions(modifiedPositions);
  r = calculator_.calculate("");
  if (!r.get<Property::SuccessfulCalculation>())
    throw std::runtime_error("Energy calculation in numerical Hessian evaluation failed.");
  auto Em = r.get<Property::Energy>();

  modifiedPositions.row(atomI)(compI) = referencePositions.row(atomI)(compI) + delta;
  calculator_.modifyPositions(std::move(modifiedPositions));
  r = calculator_.calculate("");
  if (!r.get<Property::SuccessfulCalculation>())
    throw std::runtime_error("Energy calculation in numerical Hessian evaluation failed.");
  auto Ep = r.get<Property::Energy>();

  return (Ep - 2 * E + Em) / (delta * delta);
}

double NumericalHessianCalculator::hessianElementDifferentFromEnergy(int i, int j, const PositionCollection& referencePositions,
                                                                     double delta) {
  auto modifiedPositions = referencePositions;
  auto D = delta / 2;

  auto atomI = i / 3;
  auto compI = i % 3;
  auto atomJ = j / 3;
  auto compJ = j % 3;

  modifiedPositions.row(atomI)[compI] = referencePositions.row(atomI)(compI) + D;
  modifiedPositions.row(atomJ)(compJ) = referencePositions.row(atomJ)(compJ) + D;
  calculator_.modifyPositions(modifiedPositions);
  auto r = calculator_.calculate("");
  if (!r.get<Property::SuccessfulCalculation>())
    throw std::runtime_error("Energy calculation in numerical Hessian evaluation failed.");
  auto Epp = r.get<Property::Energy>();

  modifiedPositions.row(atomI)(compI) = referencePositions.row(atomI)(compI) - D;
  modifiedPositions.row(atomJ)(compJ) = referencePositions.row(atomJ)(compJ) + D;
  calculator_.modifyPositions(modifiedPositions);
  r = calculator_.calculate("");
  if (!r.get<Property::SuccessfulCalculation>())
    throw std::runtime_error("Energy calculation in numerical Hessian evaluation failed.");
  auto Emp = r.get<Property::Energy>();

  modifiedPositions.row(atomI)(compI) = referencePositions.row(atomI)(compI) + D;
  modifiedPositions.row(atomJ)(compJ) = referencePositions.row(atomJ)(compJ) - D;
  calculator_.modifyPositions(modifiedPositions);
  r = calculator_.calculate("");
  if (!r.get<Property::SuccessfulCalculation>())
    throw std::runtime_error("Energy calculation in numerical Hessian evaluation failed.");
  auto Epm = r.get<Property::Energy>();

  modifiedPositions.row(atomI)(compI) = referencePositions.row(atomI)(compI) - D;
  modifiedPositions.row(atomJ)(compJ) = referencePositions.row(atomJ)(compJ) - D;
  calculator_.modifyPositions(modifiedPositions);
  r = calculator_.calculate("");
  if (!r.get<Property::SuccessfulCalculation>())
    throw std::runtime_error("Energy calculation in numerical Hessian evaluation failed.");
  auto Emm = r.get<Property::Energy>();

  return (Epp - Epm - Emp + Emm) / (delta * delta);
}

Results NumericalHessianCalculator::calculateFromGradientDifferences(double delta) {
  PropertyList p;
  if (dipoleGradient_) {
    p = Property::Gradients | Property::Dipole;
  }
  else {
    p = Property::Gradients;
  }

  Results results;
  calculator_.setRequiredProperties(p);
  auto positions = calculator_.getPositions();
  auto state = calculator_.getState();
  auto nDimensions = positions.size();
  HessianMatrix H = HessianMatrix::Zero(nDimensions, nDimensions);
  DipoleGradient dG = DipoleGradient::Zero(nDimensions, 3);

  std::exception_ptr exception = nullptr;
#pragma omp parallel shared(exception)
  {
    std::shared_ptr<Core::Calculator> localCalculator;
#pragma omp critical(clone)
    { localCalculator = calculator_.clone(); }
    localCalculator->setRequiredProperties(p);

#pragma omp for reduction(+ : dG)
    for (int i = 0; i < nDimensions; ++i) {
      if (!exception) {
        try {
          Eigen::VectorXd gradientDifference = addGradientContribution(dG, i, positions, delta, *localCalculator, state);
          H.col(i) = gradientDifference;
        }
        catch (const std::exception& e) {
#pragma omp critical
          { exception = std::current_exception(); }
        }
      }
    }
  }
  if (exception)
    std::rethrow_exception(exception);
  // Set a compound result with hessian matrix and, if necessary, dipole gradient.
  results.set<Property::Hessian>(0.5 * (H + H.transpose()));
  if (dipoleGradient_) {
    results.set<Property::DipoleGradient>(std::move(dG));
  }

  return results;
}

Eigen::VectorXd NumericalHessianCalculator::addGradientContribution(DipoleGradient& dipoleDiff, int i,
                                                                    const Utils::PositionCollection& referencePositions,
                                                                    double delta, Core::Calculator& calculator,
                                                                    std::shared_ptr<Core::State> state) {
  calculator.loadState(std::move(state));
  auto modifiedPositions = referencePositions;
  auto D = delta / 2;

  auto atomI = i / 3;
  auto compI = i % 3;

  modifiedPositions.row(atomI)(compI) = referencePositions.row(atomI)(compI) + D;
  calculator.modifyPositions(modifiedPositions);
  const auto& rPlus = calculator.calculate("");
  if (!rPlus.get<Property::SuccessfulCalculation>())
    throw std::runtime_error("Gradient calculation in semi-numerical Hessian evaluation failed.");
  auto gradDifference = rPlus.get<Property::Gradients>();
  if (dipoleGradient_) { // plus delta
    dipoleDiff.row(i) = rPlus.get<Property::Dipole>();
  }

  modifiedPositions.row(atomI)(compI) = referencePositions.row(atomI)(compI) - D;
  calculator.modifyPositions(std::move(modifiedPositions));
  const auto& rMinus = calculator.calculate("");
  gradDifference -= rMinus.get<Property::Gradients>();
  if (!rMinus.get<Property::SuccessfulCalculation>())
    throw std::runtime_error("Gradient calculation in semi-numerical Hessian evaluation failed.");
  if (dipoleGradient_) { // minus delta
    dipoleDiff.row(i) -= rMinus.get<Property::Dipole>();
    dipoleDiff.row(i) /= delta;
  }

  // Linearize the GradientCollection to be ready to be summed up in the hessian
  Eigen::RowVectorXd gradientDiff = Eigen::Map<Eigen::RowVectorXd>(gradDifference.data(), gradDifference.size());
  return gradientDiff / delta;
}

Results NumericalHessianCalculator::calculateFromGradientDifferences(double delta, const std::vector<int>& atomsToUpdate) {
  PropertyList p;
  if (dipoleGradient_) {
    p = Property::Gradients | Property::Dipole;
  }
  else {
    p = Property::Gradients;
  }
  calculator_.setRequiredProperties(p);
  auto state = calculator_.getState();
  const PositionCollection& positions = calculator_.getPositions();

  auto nDimensions = positions.size();
  HessianMatrix hessian = HessianMatrix::Zero(nDimensions, nDimensions);
  DipoleGradient dipoleGradient = DipoleGradient::Zero(nDimensions, 3);

  std::exception_ptr exception = nullptr;
#pragma omp parallel shared(exception)
  {
    std::shared_ptr<Core::Calculator> localCalculator;
#pragma omp critical(clone)
    { localCalculator = calculator_.clone(); }
    localCalculator->setRequiredProperties(p);
#pragma omp for reduction(+ : dipoleGradient)
    for (int i = 0; i < atomsToUpdate.size(); ++i) {
      for (int dimension = 0; dimension < 3; ++dimension) {
        if (!exception) {
          try {
            hessian.col(atomsToUpdate[i] * 3 + dimension) = addGradientContribution(
                dipoleGradient, atomsToUpdate[i] * 3 + dimension, positions, delta, *localCalculator, state);
          }
          catch (const std::exception& e) {
#pragma omp critical
            { exception = std::current_exception(); }
          }
        }
      }
    }
  }
  if (exception)
    std::rethrow_exception(exception);

  hessian += hessian.transpose().eval();
  for (int i : atomsToUpdate) {
    for (int j : atomsToUpdate) {
      hessian.block<3, 3>(i * 3, j * 3) /= 2.;
    }
  }
  Results results;
  results.set<Property::Hessian>(hessian);
  if (dipoleGradient_) {
    results.set<Property::DipoleGradient>(std::move(dipoleGradient));
  }
  return results;
}

} // namespace Utils
} // namespace Scine
