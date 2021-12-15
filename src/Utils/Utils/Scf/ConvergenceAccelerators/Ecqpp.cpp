/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Ecqpp.h"
#include <Eigen/QR>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>

namespace Scine {
namespace Utils {

Ecqpp::Ecqpp(const Eigen::MatrixXd& B, const Eigen::VectorXd& E)
  : B_(B), E_(E), dimension_(static_cast<unsigned>(E.size())) {
  assert(B.rows() == dimension_ && B.cols() == dimension_ && "B matrix and E vector do not have compatible dimensions.");
}

Eigen::VectorXd Ecqpp::calculateOptimalCoefficients() {
  solveAllConstrainedProblems();
  return solution_;
}

void Ecqpp::solveAllConstrainedProblems() {
  bestSolutionEnergy_ = std::numeric_limits<double>::max();
  for (unsigned i = 0; i < dimension_; ++i) {
    solveAllConstrainedProblemsForNumberZeros(i);
  }
}

void Ecqpp::solveAllConstrainedProblemsForNumberZeros(unsigned int numberZeros) {
  std::vector<bool> indexConsidered(dimension_, true);
  for (unsigned i = 0; i < numberZeros; ++i) {
    indexConsidered[i] = false;
  }

  do {
    generatePreviousIndexesVector(indexConsidered, numberZeros);
    generateReducedObjects();
    solveConstrainedProblem();
    if (solutionIsValid()) {
      addSolution();
    }
  } while (std::next_permutation(indexConsidered.begin(), indexConsidered.end()));
}

void Ecqpp::generatePreviousIndexesVector(const std::vector<bool>& consideredIndexes, unsigned numberZeros) {
  unsigned writeIndex = 0;
  unsigned newMatrixDimension = dimension_ - numberZeros;
  previousIndexes_.resize(newMatrixDimension);
  for (unsigned readIndex = 0; readIndex < dimension_; ++readIndex) {
    if (consideredIndexes[readIndex]) {
      previousIndexes_[writeIndex] = readIndex;
      ++writeIndex;
    }
  }
}

void Ecqpp::generateReducedObjects() {
  auto reducedDimension = static_cast<unsigned>(previousIndexes_.size());

  reducedB_.resize(reducedDimension, reducedDimension);
  reducedE_.resize(reducedDimension);

  for (unsigned i = 0; i < reducedDimension; ++i) {
    reducedE_(i) = E_(previousIndexes_[i]);
    for (unsigned j = 0; j < reducedDimension; ++j) {
      reducedB_(i, j) = B_(previousIndexes_[i], previousIndexes_[j]);
    }
  }
}

void Ecqpp::solveConstrainedProblem() {
  auto rDim = reducedB_.cols();

  Eigen::MatrixXd M(rDim + 1, rDim + 1);
  Eigen::VectorXd b(rDim + 1);

  M.block(0, 0, rDim, rDim) = reducedB_;
  M.block(rDim, 0, 1, rDim).setOnes();
  M.block(0, rDim, rDim, 1).setOnes();
  M(rDim, rDim) = 0;
  b.head(rDim) = reducedE_;
  b(rDim) = 1;

  Eigen::VectorXd solution = M.colPivHouseholderQr().solve(b);
  reducedSolution_ = solution.head(rDim);
}

bool Ecqpp::solutionIsValid() const {
  return (reducedSolution_.array() >= 0).all();
}

void Ecqpp::addSolution() {
  generateSolutionFromReducedSolution();
  setBestSolutionIfHasLowerEnergy();
}

void Ecqpp::generateSolutionFromReducedSolution() {
  fullSolution_ = Eigen::VectorXd::Zero(dimension_);
  for (unsigned i = 0; i < previousIndexes_.size(); ++i) {
    fullSolution_[previousIndexes_[i]] = reducedSolution_[i];
  }
}

void Ecqpp::setBestSolutionIfHasLowerEnergy() {
  double currentEnergy = E_.dot(fullSolution_) - 0.5 * fullSolution_.transpose() * B_ * fullSolution_;
  if (currentEnergy < bestSolutionEnergy_) {
    solution_ = fullSolution_;
    bestSolutionEnergy_ = currentEnergy;
  }
}
} // namespace Utils
} // namespace Scine
