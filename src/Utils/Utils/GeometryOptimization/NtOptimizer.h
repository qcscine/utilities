/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_NTOPTIMIZER_H_
#define UTILS_NTOPTIMIZER_H_

#include "CoordinateSystem.h"
#include "Utils/CalculatorBasics/PropertyList.h"
#include "Utils/CalculatorBasics/Results.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Geometry/ElementInfo.h"
#include "Utils/Geometry/GeometryUtilities.h"
#include "Utils/Geometry/InternalCoordinates.h"
#include "Utils/Optimizer/Optimizer.h"
#include "Utils/UniversalSettings/OptimizationSettingsNames.h"
#include <Core/Interfaces/Calculator.h>
#include <Eigen/Core>

namespace Scine {
namespace Utils {

class NtOptimizerSettings;

/**
 * @brief A version of the GeometryOptimizer that optimizes the underlying structure while applying an
 *        additional force in order to climb towards transition states.
 *
 * This version of the Newton Trajectory (NT) optimizer is inspired by the comment made in:
 * J. Comput. Chem. 2019, 9999, 1 [DOI: 10.1002/jcc.26115] ,
 * https://onlinelibrary.wiley.com/doi/full/10.1002/jcc.26115 ,
 *
 * Which references:
 * J. Comput. Chem. 2016, 37, 2467 [DOI: 10.1002/jcc.24470]
 * https://onlinelibrary.wiley.com/doi/10.1002/jcc.24470 ,
 * and:
 * Theor. Chem. Acc. 2016, 135, 113. [DOI: 10.1007/s00214-016-1880-2]
 * https://doi.org/10.1007/s00214-016-1880-2
 *
 * which has also been used as a reference.
 * However, none of the papers describe the algorithm used here in its exact form.
 */
class NtOptimizer : public Optimizer {
 public:
  /**
   * @brief Construct a new NtOptimizer object.
   * @param calculator The calculator to be used for the underlying single point/gradient calculations.
   */
  NtOptimizer(Core::Calculator& calculator) : _calculator(calculator){};
  /**
   * @brief See NtOptimizerBase::optimize().
   *
   * @param atoms The AtomCollection (Geometry) to be optimized.
   * @param log The logger to which eventual output is written.
   * @return int  The final number of optimization cycles carried out.
   */
  virtual int optimize(AtomCollection& atoms, Core::Log& log) final;
  /**
   * @brief Function to apply the given settings to underlying classes.
   * @param settings The new settings.
   */
  virtual void setSettings(const Settings& settings) final;
  void applySettings(const Settings& settings) final;
  /**
   * @brief Get the public settings as a Utils::Settings object.
   * @return Settings A settings object with the current settings.
   */
  virtual Settings getSettings() const final;
  // @brief The list indices of atoms to be forced onto or away from those in the RHS list.
  std::vector<int> lhsList = {};
  // @brief The list indices of atoms to be forced onto or away from those in the LHS list.
  std::vector<int> rhsList = {};
  /**
   * @brief The norm of the summed additional forces acting on all listed atoms.
   *
   * This value tailors the speed at which the affected atoms move in the optimization.
   * It is crucial not to move faster than the chosen optimizer can minimize all other unaffected atoms to
   * follow the trajectory.
   */
  double totalForceNorm = 0.1;
  /**
   * @brief Switch for the additional force to be attractive or repulsive.
   *
   * If true, the atoms in the lhsList and rhsList will be forced towards one another,
   * if false, the atoms will be repelled from one another.
   */
  bool attractive = true;
  /**
   * @brief Set the coordinate system in which the optimization shall be performed
   *
   * The optimization can be carried out in the internal coordinate space or with removed translations and rotations
   * possibly accelerating convergence.
   */
  CoordinateSystem coordinateSystem = CoordinateSystem::CartesianWithoutRotTrans;
  /**
   * @brief Indices of atoms that are constrained.
   */
  std::vector<int> fixedAtoms;
  /**
   * @brief Side that is moved towards/away
   */
  std::string movableSide = "both";
  /**
   * @brief If true, uses a BFGS/GDIIS in between forced steps to run some
   *        constrained geometry optimizations.
   */
  bool useMicroCycles = true;
  /**
   * @brief If true, uses `numberOfMicroCycles`, if false allows for more micro
   *        cycles as the number of NT steps grow (1 more per NT cycle).
   */
  bool fixedNumberOfMicroCycles = true;
  /**
   * @brief The fixed number of micro cycles.
   */
  int numberOfMicroCycles = 10;
  /**
   * @brief Number of passes through a Savitzky-Golay filter before analyzing
   *        the reaction curve.
   */
  int filterPasses = 10;
  // @brief possible options for extraction
  const std::vector<std::string> possibleExtractionOptions = {
      std::string(SettingsNames::Optimizations::Nt::extractHighest),
      std::string(SettingsNames::Optimizations::Nt::extractFirst),
  };
  /**
   * @brief Criterion to extract a TS guess from the trajectory
   */
  std::string extractionCriterion = possibleExtractionOptions.front();
  /// @brief The special convergence settings for this optimizer.
  struct NtConvergenceStub {
    /// @brief The maximum number of iterations.
    unsigned int maxIter = 500;
    /**
     * @brief The minimum distance, given as multiple of covalent radii sums,
     *        between all atoms in the constrained lists in case of an attractive
     *        run.
     */
    double attractiveStop = 0.9;
    /**
     * @brief The distance, given as multiple of covalent radii sums, all atoms
     *        in the constrained lists have to be apart in case of a repulsive
     *        run.
     */
    double repulsiveStop = 4.0;
  };
  /// @brief The convergence control for this optimizer.
  NtConvergenceStub check;
  /// @brief The special settings for the optimizer of the macro iterations.
  struct NtOptimizerStub {
    /// @brief The steepest descent factor.
    double factor = 1.0;
  };
  /// @brief The optimizer for the macro iterations.
  NtOptimizerStub optimizer;
  /**
   * @brief @see Scine::Utils::Optimizer::addSettingsDescriptors
   * @param collection The collection to add the descriptors to.
   */
  void addSettingsDescriptors(UniversalSettings::DescriptorCollection& collection) const final;

 private:
  /**
   * @brief Checks the lhs and rhs list for valid indices
   *
   * @throw std::logic_error if invalid index
   */
  void sanityCheck(const AtomCollection& atoms) const;
  /**
   * @brief Calculates vector from geometric centers of lhs and rhs.
   *
   * @param positions The current positions.
   * @return Displacement The vector from rhs to lhs
   */
  Displacement centerToCenterVector(const PositionCollection& positions) const;
  /**
   * @brief The function evaluating all artificial forces for the given set of atoms.
   *
   * @param atoms The current atoms.
   * @param energy The current energy in hartree.
   * @param gradients The gradient to be updated (in a.u).
   * @param addForce Add external force.
   */
  void updateGradients(const AtomCollection& atoms, const double& energy, GradientCollection& gradients,
                       bool addForce = false) const;

  /**
   * @brief The function determining whether the optimization is converged.
   *
   * @param atoms The current atoms.
   * @return bool whether it is converged
   */
  bool convergedOptimization(const AtomCollection& atoms) const;
  /**
   * @brief Extracts the transition state (TS) guess from the generated trajectory.
   * @return PositionCollection The transition state guess.
   */
  PositionCollection extractTsGuess() const;
  /**
   * @brief The function performing the SD step on the positions.
   *
   * @param coordinates The current positions, which are changed in place.
   * @param atoms The current atoms to perform coordinate transformations.
   * @param gradients The current gradients for the SD step.
   */
  void updateCoordinates(PositionCollection& coordinates, const AtomCollection& atoms, const GradientCollection& gradients) const;
  // @brief values of macro cycles
  std::vector<double> _values;
  // @brief trajectory of macro cycles
  std::vector<PositionCollection> _trajectory;
  Core::Calculator& _calculator;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_NTOPTIMIZER_H_
