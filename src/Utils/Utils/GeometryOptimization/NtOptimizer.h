/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_NTOPTIMIZER_H_
#define UTILS_NTOPTIMIZER_H_

#include "Utils/CalculatorBasics/PropertyList.h"
#include "Utils/CalculatorBasics/Results.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Geometry/ElementInfo.h"
#include "Utils/Geometry/GeometryUtilities.h"
#include "Utils/Geometry/InternalCoordinates.h"
#include "Utils/Optimizer/Optimizer.h"
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
  // Definition of Utils::Settings keys
  static constexpr const char* ntRHSListKey = "nt_rhs_list";
  static constexpr const char* ntLHSListKey = "nt_lhs_list";
  static constexpr const char* ntAttractiveKey = "nt_attractive";
  static constexpr const char* ntTotalForceNormKey = "nt_total_force_norm";
  static constexpr const char* ntTransfromCoordinatesKey = "nt_transform_coordinates";
  static constexpr const char* ntMaxIterKey = "convergence_max_iterations";
  static constexpr const char* ntRepulsiveStopKey = "convergence_repulsive_stop";
  static constexpr const char* ntAttractiveStopKey = "convergence_attractive_stop";
  static constexpr const char* ntSdFactorKey = "sd_factor";
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
  virtual void applySettings(const Settings& settings) final;
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
  double totalForceNorm = 0.01;
  /**
   * @brief Switch for the additional force to be attractive or repulsive.
   *
   * If true, the atoms in the lhsList and rhsList will be forced towards one another,
   * if false, the atoms will be repelled from one another.
   */
  bool attractive = true;
  /**
   * @brief Switch to transform the coordinates from Cartesian into an internal space.
   *
   * The optimization will be carried out in the internal coordinate space possibly
   * accelerating convergence.
   */
  bool transformCoordinates = true;
  /// @brief The special convergence settings for this optimizer.
  struct NtConvergenceStub {
    /// @brief The maximum number of iterations.
    unsigned int maxIter = 200;
    /**
     * @brief The minimum distance, given as multiple of covalent radii sums,
     *        between all atoms in the contrained lists in case of an attractive
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
  /// @brief The special convergence settings for this optimizer.
  struct NtOptimizerStub {
    /// @brief The maximum number of iterations.
    double factor = 1.0;
  };
  /// @brief The convergence control for this optimizer.
  NtOptimizerStub optimizer;
  /**
   * @brief @see Scine::Utils::Optimizer::addSettingsDescriptors
   * @param collection The collection to add the descriptors to.
   */
  virtual void addSettingsDescriptors(UniversalSettings::DescriptorCollection& collection) const final;

 private:
  /**
   * @brief The function evaluating all artificial forces for the given set of atoms.
   *
   * @param atoms The current atom positions.
   * @param energy The current energy in hartree.
   * @param gradients The gradient to be updated (in a.u).
   * @return double The norm of the real gradient projected onto the applied forces direction.
   */
  double updateGradients(const AtomCollection& atoms, const double& energy, GradientCollection& gradients) const;
  /**
   * @brief Extracts the transition state (TS) guess from the generated trajectory.
   * @return PositionCollection The transition state guess.
   */
  PositionCollection extractTsGuess() const;
  std::vector<double> _norms;
  std::vector<PositionCollection> _trajectory;
  Core::Calculator& _calculator;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_NTOPTIMIZER_H_
