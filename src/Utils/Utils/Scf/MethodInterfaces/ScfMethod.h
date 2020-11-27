/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_SCFMETHOD_H
#define UTILS_SCFMETHOD_H

#include "ConvergenceChecker.h"
#include "LcaoMethod.h"
#include "ScfConvergenceAccelerator.h"
#include <Utils/Math/DerivOrderEnum.h>
#include <map>

namespace Scine {
namespace Utils {

class ScfModifier;
class DensityMatrixGuessCalculator;
enum class scf_mixer_t;

/**
 * @class ScfMethod @file ScfMethod.h
 * @brief Base class for self-consistent-field methods.
 */
class ScfMethod : public LcaoMethod {
 public:
  explicit ScfMethod(bool unrestrictedCalculationPossible, Utils::DerivativeOrder maximalOrder, bool orthogonalBasisSet = false);
  ~ScfMethod() override;

  void initialize() override;

  void calculate(Utils::Derivative d, Core::Log& log) override;

  /*! Adds a modifier to the SCF method.
      First checks whether the ScfModifier already exists.
      Will call the functions setMethod(...) and initialize of modifier.
      @param priority Number from 0 (early) to 10 (late) to set when the modifier will be called compared to other
     ones.*/
  void addModifier(std::shared_ptr<ScfModifier> modifier, int priority = 5);
  /*! Removes the SCF modifier. */
  void removeModifier(const std::shared_ptr<ScfModifier>& modifier);

  /*! Performs an SCF calculation (NB: will stop after max number of cycles). */
  void convergedCalculation(Core::Log& log, Utils::Derivative d = Utils::Derivative::First);
  void performIteration(Core::Log& log, Utils::Derivative d = Utils::Derivative::First);

  bool hasConverged() const {
    return converged;
  }
  /*! Returns the number of iterations that were performed for the energy calculation. */
  int getNumberIterations() const {
    return performedIterations_;
  }

  void setScfMixer(scf_mixer_t mixer);
  scf_mixer_t getScfMixer() const;
  /*! Set the convergence threshold for a SCF calculation (RMSD of density matrix) */
  void setConvergenceCriteria(double c) {
    convergenceChecker_.setThreshold(c);
  }
  /*! Get the currently applied threshold for SCF calculation. */
  double getConvergenceThreshold() const {
    return convergenceChecker_.getThreshold();
  }
  void setMaxIterations(int max) {
    maxIterations = max;
  }

  /*! Reset the convergence check; is useful if the starting density matrix in the SCF cycle is not the last one of the
   * previous cycle. */
  void resetConvergenceCheck();
  int getMaxNumberIterations() const {
    return maxIterations;
  }

  /*! Reset the density matrix. */
  void reinitializeDensityMatrixGuess();
  /*! Returns the density matrix guess for the current system */
  DensityMatrix getDensityMatrixGuess() const;

  /*! @name Steps of SCF calculations
   * @{
   */
  void calculateDensityDependentQuantities(Utils::Derivative d = Utils::Derivative::First);
  void finalizeCalculation(Utils::Derivative d = Utils::Derivative::First);
  /*! @} */

  /*! Perform a calculation for the current density matrix without optimizing it. */
  void evaluateDensity(Utils::Derivative derivativeOrder);

 protected:
  std::shared_ptr<DensityMatrixGuessCalculator> densityMatrixGuess_;
  bool converged;
  int performedIterations_; // Tells how many iterations were performed
  int maxIterations;

 private:
  void solveEigenValueProblem(Core::Log& log);
  void onConvergedCalculationStarts();
  using ModifierContainer = std::multimap<int, std::shared_ptr<ScfModifier>>;

  ModifierContainer modifiers;
  ConvergenceChecker convergenceChecker_;
  ScfConvergenceAccelerator convergenceAccelerator_;
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_SCFMETHOD_H
