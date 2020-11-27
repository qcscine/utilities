/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_AFIROPTIMIZERBASE_H_
#define UTILS_AFIROPTIMIZERBASE_H_

#include "Utils/Typenames.h"
#include <Core/Log.h>
#include <Eigen/Core>
#include <vector>

namespace Scine {
namespace Utils {

class AtomCollection;
class Settings;

/**
 * @brief The base class for all AFIR optimizers.
 *
 * The main purpose of this base class is to hide the template parameter(s)
 * of the derived class(es).
 *
 * This optimizer is a version of the GeometryOptimizer that optimizes the underlying structure
 * while applying an additional artificial force for the induction of reactions.
 *
 * The artificial force induced reaction (AFIR) optimizer is based on the works by Maeda et. al.
 * The implementation was done following this reference:
 * J Comput Chem., 2018, 39, 233 [DOI: 10.1002/jcc.25106]
 * https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5765425/ (accessed 05.04.2019)
 *
 * The original references are:
 * J. Chem. Phys., 2010, 132, 241102
 * J. Chem. Theory Comput., 2011, 7, 2335
 */
class AfirOptimizerBase {
 public:
  // Definition of Utils::Settings keys
  static constexpr const char* afirRHSListKey = "afir_rhs_list";
  static constexpr const char* afirLHSListKey = "afir_lhs_list";
  static constexpr const char* afirWeakForcesKey = "afir_weak_forces";
  static constexpr const char* afirAttractiveKey = "afir_attractive";
  static constexpr const char* afirEnergyAllowanceKey = "afir_energy_allowance";
  static constexpr const char* afirPhaseInKey = "afir_phase_in";
  static constexpr const char* afirTransfromCoordinatesKey = "afir_transform_coordinates";

  /// @brief Default constructor.
  AfirOptimizerBase() = default;
  /// @brief Virtual default destructor.
  virtual ~AfirOptimizerBase() = default;
  /**
   * @brief The main functionality of the geometry optimizer.
   *
   * This function wraps the optimize functions of the underlying optimizer.
   *
   * @param atoms The AtomCollection (Geometry) to be optimized.
   * @param log The logger to which eventual output is written.
   * @return int  The final number of optimization cycles carried out.
   */
  virtual int optimize(AtomCollection& atoms, Core::Log& log) = 0;
  /**
   * @brief Function to apply the given settings to underlying classes.
   * @param settings The new settings.
   */
  virtual void setSettings(const Settings& settings) = 0;
  /**
   * @brief Get the public settings as a Utils::Settings object.
   * @return Settings A settings object with the current settings.
   */
  virtual Settings getSettings() const = 0;
  /**
   * @brief Add an observer function that will be triggered in each iteration.
   *
   * @param function A function to be executed in every loop of the optimization.
   *                 The function will have access to the current cycle count,
   *                 the current value and to a const reference of the current
   *                 parameters.
   */
  virtual void addObserver(std::function<void(const int&, const double&, const Eigen::VectorXd&)> function) = 0;
  /**
   * @brief Clear all existing observer functions.
   *
   * For optimization problems with very fast evaluations of the underlying function
   * the removal of all observers can increase performance as the observers are given
   * as std::functions and can not be added via templates.
   */
  virtual void clearObservers() = 0;
  /// @brief The list of atoms indices of atoms to be artificially forced onto or away from those in the RHS list.
  std::vector<int> lhsList = {};
  /// @brief The list of atoms indices of atoms to be artificially forced onto or away from those in the LHS list.
  std::vector<int> rhsList = {};
  /**
   * @brief The maximum amount of energy to be added by the artifical force, in kJ/mol.
   *
   * This values gives the maximum height of an energy barrier the artificial force can force the molecule
   * across.
   */
  double energyAllowance = 1000.0;
  /**
   * @brief Switch for the artificial force to be attractive or repulsive.
   *
   * If true, the atoms in the lhsList and rhsList will be forced towards one another,
   * if false, the atoms will be repelled from one another.
   */
  bool attractive = true;
  /**
   * @brief Switch for an additional weak attractive force applied to all atom pairs.
   *
   * If true a small attractive force is added across the entire set of atoms.
   * The energy allowance for this force is calculated as \f$ 10/(N_{\text{atoms}}-1) \f$ kJ/mol.
   */
  bool weak = false;
  /**
   * @brief The number of steps over which the full attractive force is slowly applied.
   *
   * The phasing in of the artificial forces happens linearly.
   * Note that only the force (gradient contribution) is phased in, the energy contribution
   * is always added fully.
   */
  int phaseIn = 100;
  /**
   * @brief Switch to transform the coordinates from Cartesian into an internal space.
   *
   * The optimization will be carried out in the internal coordinate space possibly
   * accellerating convergence.
   */
  bool transformCoordinates = true;
  /**
   * @brief The function evaluating all artificial forces for the given set of atoms.
   *
   * @param atoms The atoms.
   * @param energy The added artificial energy in hartree/a.u.
   * @param gradients The added artificial gradient contributions in a.u.
   */
  void evaluateArtificialForces(const AtomCollection& atoms, double& energy, GradientCollection& gradients) const;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_AFIROPTIMIZERBASE_H_
