/**
 * @file NumericalHessianCalculator.cpp
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "NumericalHessianCalculator.h"
#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/CalculatorBasics/StatesHandler.h>

namespace Scine {
namespace Utils {

#pragma omp declare reduction(+ : Eigen::MatrixXd : omp_out = omp_out + omp_in) initializer(omp_priv = omp_orig)
#pragma omp declare reduction(+ : DipoleGradient : omp_out = omp_out + omp_in) initializer(omp_priv = omp_orig)

NumericalHessianCalculator::NumericalHessianCalculator(Core::Calculator& calculator) : calculator_(calculator) {
}

void NumericalHessianCalculator::requiredDipoleGradient(bool dipoleGradient) {
  dipoleGradient_ = dipoleGradient;
}

Results NumericalHessianCalculator::calculate(double delta) {
  return calculateFromGradientDifferences(delta);
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
  auto E = r.getEnergy();

  modifiedPositions.row(atomI)(compI) = referencePositions.row(atomI)(compI) - delta;
  calculator_.modifyPositions(modifiedPositions);
  r = calculator_.calculate("");
  auto Em = r.getEnergy();

  modifiedPositions.row(atomI)(compI) = referencePositions.row(atomI)(compI) + delta;
  calculator_.modifyPositions(std::move(modifiedPositions));
  r = calculator_.calculate("");
  auto Ep = r.getEnergy();

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
  auto Epp = r.getEnergy();

  modifiedPositions.row(atomI)(compI) = referencePositions.row(atomI)(compI) - D;
  modifiedPositions.row(atomJ)(compJ) = referencePositions.row(atomJ)(compJ) + D;
  calculator_.modifyPositions(modifiedPositions);
  r = calculator_.calculate("");
  auto Emp = r.getEnergy();

  modifiedPositions.row(atomI)(compI) = referencePositions.row(atomI)(compI) + D;
  modifiedPositions.row(atomJ)(compJ) = referencePositions.row(atomJ)(compJ) - D;
  calculator_.modifyPositions(modifiedPositions);
  r = calculator_.calculate("");
  auto Epm = r.getEnergy();

  modifiedPositions.row(atomI)(compI) = referencePositions.row(atomI)(compI) - D;
  modifiedPositions.row(atomJ)(compJ) = referencePositions.row(atomJ)(compJ) - D;
  calculator_.modifyPositions(modifiedPositions);
  r = calculator_.calculate("");
  auto Emm = r.getEnergy();

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
  auto state = calculator_.statesHandler().getCurrentState(StateSize::minimal);
  auto nDimensions = positions.size();
  HessianMatrix H = HessianMatrix::Zero(nDimensions, nDimensions);
  DipoleGradient dG;
  dG.resize(nDimensions, 3);
  std::vector<Eigen::RowVectorXd> gradientDifferences(nDimensions, Eigen::RowVectorXd(nDimensions));
  std::vector<Dipole> dipoleDifferences;
  if (dipoleGradient_) {
    dipoleDifferences.resize(nDimensions);
  }

#pragma omp parallel
  {
    std::shared_ptr<Core::Calculator> localCalculator;
#pragma omp critical(clone)
    { localCalculator = calculator_.clone(); }
    localCalculator->setRequiredProperties(p);

#pragma omp for reduction(+ : dG)
    for (int i = 0; i < nDimensions; ++i) {
      Eigen::VectorXd gradientDifference = addGradientContribution(dG, i, positions, delta, localCalculator, state);
      H.col(i) = gradientDifference;
    }
  }
  // Set a compound result with hessian matrix and, if necessary, dipole gradient.
  results.setHessian(H);
  if (dipoleGradient_) {
    dG /= delta;
    results.setDipoleGradient(std::move(dG));
  }

  return results;
}

Eigen::VectorXd NumericalHessianCalculator::addGradientContribution(DipoleGradient& dipoleDiff, int i,
                                                                    const Utils::PositionCollection& referencePositions,
                                                                    double delta, std::shared_ptr<Core::Calculator> calculator,
                                                                    std::shared_ptr<State> state) {
  calculator->statesHandler().load(std::move(state));
  auto modifiedPositions = referencePositions;
  auto D = delta / 2;

  auto atomI = i / 3;
  auto compI = i % 3;

  modifiedPositions.row(atomI)(compI) = referencePositions.row(atomI)(compI) + D;
  calculator->modifyPositions(modifiedPositions);
  auto rPlus = calculator->calculate("");
  auto gradDifference = rPlus.takeGradients();
  if (dipoleGradient_) { // plus delta
    dipoleDiff.row(atomI) = rPlus.getDipole();
  }

  modifiedPositions.row(atomI)(compI) = referencePositions.row(atomI)(compI) - D;
  calculator->modifyPositions(std::move(modifiedPositions));
  auto rMinus = calculator->calculate("");
  gradDifference -= rMinus.takeGradients();
  if (dipoleGradient_) { // minus delta
    dipoleDiff.row(atomI) -= rMinus.getDipole();
  }

  // Linearize the GradientCollection to be ready to be summed up in the hessian
  auto gradientDiff = Eigen::Map<Eigen::RowVectorXd>(gradDifference.data(), gradDifference.size());
  return gradientDiff / delta;
}

} // namespace Utils
} // namespace Scine
