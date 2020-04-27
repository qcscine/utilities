/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_EDIIS_ECQPP_H
#define UTILS_EDIIS_ECQPP_H

#include <Eigen/Core>
#include <vector>

namespace Scine {
namespace Utils {

/*!
 * Equality Constrained Quadratic Programming Problem solver for EDIIS.
 * This class solves the optimization problem for an EDIIS iteration.
 * It minimizes the expression E*c - 0.5 c^T*B*c with the constraint that sum(c_i) = 1.
 * To do so, it solves the above equation for all subsets of c, rejects inadmissible solutions (c_i < 0) and keeps the
 * best one.
 */
class Ecqpp {
 public:
  Ecqpp(const Eigen::MatrixXd& B, const Eigen::VectorXd& E);

  Eigen::VectorXd calculateOptimalCoefficients();

 private:
  void solveAllConstrainedProblems();
  void solveAllConstrainedProblemsForNumberZeros(unsigned int numberZeros);
  void generatePreviousIndexesVector(const std::vector<bool>& consideredIndexes, unsigned numberZeros);
  void generateReducedObjects();
  void solveConstrainedProblem();
  bool solutionIsValid() const;
  void addSolution();
  void generateSolutionFromReducedSolution();
  void setBestSolutionIfHasLowerEnergy();

  const Eigen::MatrixXd& B_;
  const Eigen::VectorXd& E_;
  const unsigned dimension_;

  std::vector<unsigned> previousIndexes_;
  Eigen::MatrixXd reducedB_;
  Eigen::VectorXd reducedE_;
  Eigen::VectorXd reducedSolution_;
  Eigen::VectorXd fullSolution_;
  Eigen::VectorXd solution_;
  double bestSolutionEnergy_;
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_EDIIS_ECQPP_H