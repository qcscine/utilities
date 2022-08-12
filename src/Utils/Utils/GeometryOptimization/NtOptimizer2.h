/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_NTOPTIMIZER2_H_
#define UTILS_NTOPTIMIZER2_H_

#include "CoordinateSystem.h"
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

class NtOptimizer2Settings;

/**
 * @brief A version of the GeometryOptimizer that optimizes the underlying structure while applying an
 *        additional force in order to climb towards transition states.
 *
 * This version of the Newton Trajectory (NT) optimizer is a modification of the
 * first one.
 *
 * It uses the explicit information about bonds (bond orders) to define reaction
 * coordinates, it can associate atom pairs into bonds and at the same time
 * dissociate atom pairs (break bonds).
 *
 */
class NtOptimizer2 : public Optimizer {
  using ReactionMapping = std::vector<std::pair<std::vector<int>, std::vector<int>>>;

 public:
  // Definition of Utils::Settings keys
  static constexpr const char* ntAssListKey = "nt_associations";
  static constexpr const char* ntDissListKey = "nt_dissociations";
  static constexpr const char* ntTotalForceNormKey = "nt_total_force_norm";
  static constexpr const char* ntMaxIterKey = "convergence_max_iterations";
  static constexpr const char* ntRepulsiveStopKey = "convergence_repulsive_stop";
  static constexpr const char* ntAttractiveStopKey = "convergence_attractive_stop";
  static constexpr const char* ntSdFactorKey = "sd_factor";
  static constexpr const char* ntUseMicroCycles = "nt_use_micro_cycles";
  static constexpr const char* ntFixedNumberOfMicroCycles = "nt_fixed_number_of_micro_cycles";
  static constexpr const char* ntNumberOfMicroCycles = "nt_number_of_micro_cycles";
  static constexpr const char* ntFilterPasses = "nt_filter_passes";
  static constexpr const char* ntExtractionCriterion = "nt_extraction_criterion";
  static constexpr const char* ntCoordinateSystemKey = "nt_coordinate_system";
  static constexpr const char* ntFixedAtomsKey = "nt_constrained_atoms";
  static constexpr const char* ntExtractLastBeforeTarget = "last_maximum_before_first_target";
  static constexpr const char* ntExtractHighest = "highest_maximum";
  static constexpr const char* ntExtractFirst = "first_maximum";
  /**
   * @brief Construct a new NtOptimizer object.
   * @param calculator The calculator to be used for the underlying single point/gradient calculations.
   */
  NtOptimizer2(Core::Calculator& calculator) : _calculator(calculator){};
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
  /**
   * @brief The list indices of atoms to be forced onto or away from those in the RHS list.
   */
  std::vector<int> associationList = {};
  /**
   * @brief The list indices of atoms to be forced onto or away from those in the LHS list.
   */
  std::vector<int> dissociationList = {};
  /**
   * @brief Get the map of constraints to the involved atoms.
   */
  const std::vector<std::vector<int>>& getConstraintsMap();
  /**
   * @brief Get the reactive atoms list.
   */
  const std::vector<int>& getReactiveAtomsList();
  /**
   * @brief The norm of the summed additional forces acting on all listed atoms.
   *
   * This value tailors the speed at which the affected atoms move in the optimization.
   * It is crucial not to move faster than the chosen optimizer can minimize all other unaffected atoms to
   * follow the trajectory.
   */
  double totalForceNorm = 0.1;
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
   * @brief Number of passes through a Savitzky-Golay filter before analyzing the reaction curve.
   */
  int filterPasses = 10;
  // @brief possible options for extraction
  const std::vector<std::string> possibleExtractionOptions = {
      ntExtractLastBeforeTarget,
      ntExtractHighest,
      ntExtractFirst,
  };
  /**
   * @brief Criterion to extract a TS guess from the trajectory
   */
  std::string extractionCriterion = possibleExtractionOptions.front();
  /**
   * @brief The special convergence settings for this optimizer.
   */
  struct NtConvergenceStub {
    /**
     * @brief The maximum number of iterations.
     */
    unsigned int maxIter = 500;
    /**
     * @brief The minimum distance, given as multiple of covalent radii sums,
     *        between all atoms in the constrained lists in case of an attractive run.
     */
    double attractiveDistanceStop = 0.9;
    /**
     * @brief The minimum bond order necessary to signal a bond formation along an attractive reaction coordinate.
     */
    double attractiveBondOrderStop = 0.75;
    /**
     * @brief The maximum bond order necessary to still signal an existing bond along a repulsive reaction coordinate.
     */
    double repulsiveBondOrderStop = 0.15;
  };
  /**
   * @brief The convergence control for this optimizer.
   */
  NtConvergenceStub check;
  /**
   * @brief The special settings for the optimizer of the macro iterations.
   */
  struct NtOptimizerStub {
    /**
     * @brief The steepest descent factor.
     */
    double factor = 1.0;
  };
  /**
   * @brief The optimizer for the macro iterations.
   */
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
   * @brief The function evaluating all artificial forces for the given set of atoms.
   *
   * @param atoms The current atoms.
   * @param energy The current energy in hartree.
   * @param gradients The gradient to be updated (in a.u).
   * @param bondOrders The bond orders used to determine constraints
   * @param addForce Add external force.
   */
  void updateGradients(const AtomCollection& atoms, const double& energy, GradientCollection& gradients,
                       const BondOrderCollection& bondOrders, int cycle, bool addForce = false);
  /**
   * @brief Avoids gradient manipulations of the reactive atoms
   *
   * @param positions The current positions.
   * @param gradients The to be altered gradients.
   */
  void eliminateReactiveAtomsGradients(const PositionCollection& positions, GradientCollection& gradients) const;
  /**
   * @brief The function determining whether the optimization is converged.
   *
   * @param atoms The current atoms.
   * @param bondOrders The bond orders used to determine the convergence.
   * @return bool whether it is converged
   */
  bool convergedOptimization(const AtomCollection& atoms, const BondOrderCollection& bondOrders) const;
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
  /**
   * @brief Extracts the map for each associations and dissociations including eta bond considerations.
   *
   * @param bondOrders The bonds defining separate molecules.
   *
   * @return reactionMaps The reaction mapping of both associative and dissociative reaction coordinates
   */
  std::pair<ReactionMapping, ReactionMapping> inferReactions(const BondOrderCollection& bondOrders) const;
  /**
   * @brief Set the reactive atoms list by combining the sorted association and dissociation list and
   * make the combination unique.
   */
  void setReactiveAtomsList();
  /*
   * @brief Map the constraints to the involved atoms.
   *
   * @param atoms The current atoms to resize the constraints mapping list.
   */
  void setConstraintsMap(const AtomCollection& atoms);
  // @brief values of macro cycles
  std::vector<double> _values;
  // @brief trajectory of macro cycles
  std::vector<PositionCollection> _trajectory;
  // @brief A unique list of the associationList and dissociationList
  std::vector<int> _reactiveAtomsList = {};
  // @brief A mapping of atom indices to the constraints they are involved in.
  std::vector<std::vector<int>> _constraintsMap;
  Core::Calculator& _calculator;
  int _firstCoordinateReachedIndex = -1;
};

namespace NtUtils {
/**
 * @brief Gives list of list, with each sublist containing the values of indices that are part of the same molecules
 * based on the given bond orders.
 *
 * @param indices The nuclei indices that shall be sorted into molecules.
 * @param bondOrders The bonds defining separate molecules.
 *
 * @return molecules Each sublist representing a molecule.
 */
std::vector<std::vector<int>> connectedNuclei(std::vector<int> indices, const BondOrderCollection& bondOrders);
/**
 * @brief Calculates vector from geometric centers of lhsList and rhsList.
 *
 * @note opposite direction to identical method in NtOptimizer due to legacy reasons.
 *
 * @param positions The current positions.
 * @return Displacement The vector from lhs to rhs
 */
Displacement centerToCenterVector(const PositionCollection& positions, const std::vector<int>& lhsList,
                                  const std::vector<int>& rhsList);
/**
 * @brief Returns the smallest covalent radius from the given indices and atoms.
 *
 * @param atoms The atoms to pick from.
 * @param indices The list of indices of possible atoms.
 *
 * @return covRad The smallest covalent radius in Bohr
 */
double smallestCovalentRadius(const AtomCollection& atoms, const std::vector<int>& indices);

} // namespace NtUtils
} // namespace Utils
} // namespace Scine

#endif // UTILS_NTOPTIMIZER2_H_
