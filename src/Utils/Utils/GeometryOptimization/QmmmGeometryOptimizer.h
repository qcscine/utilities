/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_QMMMGEOMETRYOPTIMIZER_H
#define UTILS_QMMMGEOMETRYOPTIMIZER_H

#include "GeometryOptimizer.h"
#include "Utils/UniversalSettings/OptimizationSettingsNames.h"
#include <Core/Log.h>

namespace Scine {
namespace Utils {
/* Forward declaration */
template<class OptimizerType>
class QmmmGeometryOptimizer;

/**
 * @class QmmmGeometryOptimizerSettings QmmmGeometryOptimizer.h
 * @brief Settings for a QmmmGeometryOptimizer
 *
 * Uses template arguments in order to automatically include the
 * settings of underlying objects into the given settings.
 *
 * @tparam OptimizerType The underlying Optimizer class.
 */
template<class OptimizerType>
class QmmmGeometryOptimizerSettings : public Settings {
 public:
  /**
   * @brief Constructor.
   *
   * @param qmmmOptimizer The QM/MM geometry optimizer.
   * @param fullOptimizer The geometry optimizer for the full system.
   * @param mmOptimizer The geometry optimizer for the environment.
   */
  QmmmGeometryOptimizerSettings(const QmmmGeometryOptimizer<OptimizerType>& qmmmOptimizer,
                                const GeometryOptimizer<OptimizerType>& fullOptimizer,
                                const GeometryOptimizer<OptimizerType>& /* mmOptimizer */)
    : Settings("QmmmGeometryOptimizerSettings") {
    /*
     * The full system optimizer and the environment-only optimizers always have the same convergence criteria and
     * the same settings for the optimization algorithm (e.g., BFGS settings).
     *
     * The setting "max_iter" is overwritten by the special settings regarding the maximum number of iterations
     * that are present in the QM/MM optimizer settings. It is therefore meaningless to modify that setting
     * in the convergence check.
     */
    fullOptimizer.optimizer.addSettingsDescriptors(this->_fields);
    fullOptimizer.check.addSettingsDescriptors(this->_fields);

    UniversalSettings::OptionListDescriptor geooptCoordinateSystem("Set the coordinate system.");
    geooptCoordinateSystem.addOption("internal");
    geooptCoordinateSystem.addOption("cartesianWithoutRotTrans");
    geooptCoordinateSystem.addOption("cartesian");
    geooptCoordinateSystem.setDefaultOption(
        CoordinateSystemInterpreter::getStringFromCoordinateSystem(qmmmOptimizer.coordinateSystem));
    this->_fields.push_back(SettingsNames::Optimizations::GeometryOptimizer::coordinateSystem,
                            std::move(geooptCoordinateSystem));

    UniversalSettings::IntDescriptor qmmmOptMaxMacroiterations("The maximum number of macrocycles allowed.");
    qmmmOptMaxMacroiterations.setDefaultValue(qmmmOptimizer.maxMacrocycles);
    this->_fields.push_back(QmmmGeometryOptimizer<OptimizerType>::qmmmOptMaxMacroiterationsKey,
                            std::move(qmmmOptMaxMacroiterations));

    UniversalSettings::IntDescriptor qmmmOptMaxFullMicroiterations(
        "The maximum number of full optimization microcycles allowed per macrocycle.");
    qmmmOptMaxFullMicroiterations.setDefaultValue(qmmmOptimizer.maxFullOptMicrocycles);
    this->_fields.push_back(QmmmGeometryOptimizer<OptimizerType>::qmmmOptMaxFullMicroiterationsKey,
                            std::move(qmmmOptMaxFullMicroiterations));

    UniversalSettings::IntDescriptor qmmmOptMaxEnvMicroiteration(
        "The maximum number of MM-only optimization microcycles allowed per macrocycle.");
    qmmmOptMaxEnvMicroiteration.setDefaultValue(qmmmOptimizer.maxEnvOptMicrocycles);
    this->_fields.push_back(QmmmGeometryOptimizer<OptimizerType>::qmmmOptMaxEnvMicroiterationsKey,
                            std::move(qmmmOptMaxEnvMicroiteration));

    UniversalSettings::IntDescriptor qmmmOptEnvOptSwitchOff(
        "The number of macrocycles after which the MM-only optimization is switched off.");
    qmmmOptEnvOptSwitchOff.setDefaultValue(qmmmOptimizer.envOptSwitchOff);
    this->_fields.push_back(QmmmGeometryOptimizer<OptimizerType>::qmmmOptEnvOptSwitchOffKey, std::move(qmmmOptEnvOptSwitchOff));

    UniversalSettings::DoubleDescriptor qmmmOptBoundaryDistanceThreshold(
        "The distance threshold in Angstrom determining which atoms of the environment are also frozen during the "
        "MM-only micro-opt.");
    qmmmOptBoundaryDistanceThreshold.setDefaultValue(qmmmOptimizer.boundaryDistanceThreshold * Utils::Constants::angstrom_per_bohr);
    this->_fields.push_back(QmmmGeometryOptimizer<OptimizerType>::qmmmOptBoundaryDistanceThresholdKey,
                            std::move(qmmmOptBoundaryDistanceThreshold));

    UniversalSettings::BoolDescriptor qmmmOptEnvOptAtStart(
        "Perform the MM-only optimization of the environment at the very beginning.");
    qmmmOptEnvOptAtStart.setDefaultValue(qmmmOptimizer.envOptAtStart);
    this->_fields.push_back(QmmmGeometryOptimizer<OptimizerType>::qmmmOptEnvOptAtStartKey, std::move(qmmmOptEnvOptAtStart));

    this->resetToDefaults();
  }
};

/**
 * @class QmmmGeometryOptimizer QmmmGeometryOptimizer.h
 *
 * @brief Implementation of a QM/MM geometry optimizer that should be preferably applied for QM/MM geometry
 *        optimizations as it employs the microiteration-optimization scheme leading to a more efficient optimization
 *        with regards to the number of QM gradients calculations needed.
 *
 * Basic algorithm:
 * 1. optimize only the environment (MM region) until convergence or maximum number of
 *    iterations reached (setting: qmmm_opt_max_env_microiterations).
 * 2. optimize the full system with the QM/MM method until convergence or maximum number of
 *    iterations reached (setting: qmmm_opt_max_full_microiterations).
 * 3. repeat steps 1 and 2 (this is equal to one macrocycle) until step 2 converges. If this does not occur
 *    within the maximum number of macrocycles (setting: qmmm_opt_max_macroiterations),
 *    the whole optimization remains unconverged.
 *
 * The convergence settings and optimizer settings are the same for both geometry optimizers applied here. These
 * can be modified by the usual settings names. However, the max. iter. settings have special settings names
 * (mentioned in algorithm description above) as they have the important role of controlling the QM/MM optimization
 * algorithm. They can be set to different values for each of the different types of cycles.
 * These settings and their names are public members of this class.
 *
 * @tparam OptimizerType Expects any of the Optimizer classes. Note that some special optimizers
 *                       may not yet be supported or may need additional specialization.
 */
template<class OptimizerType>
class QmmmGeometryOptimizer {
 public:
  static constexpr const char* qmmmOptMaxMacroiterationsKey = "qmmm_opt_max_macroiterations";
  static constexpr const char* qmmmOptMaxFullMicroiterationsKey = "qmmm_opt_max_full_microiterations";
  static constexpr const char* qmmmOptMaxEnvMicroiterationsKey = "qmmm_opt_max_env_microiterations";
  static constexpr const char* qmmmOptBoundaryDistanceThresholdKey = "qmmm_opt_boundary_distance_thresh";
  static constexpr const char* qmmmOptEnvOptSwitchOffKey = "qmmm_opt_env_switch_off";
  static constexpr const char* qmmmOptEnvOptAtStartKey = "qmmm_opt_env_start";
  /**
   * @brief Constructor.
   * @param calculator The calculator to be used for the single point/gradient calculations.
   *                   The calculator has to be a QM/MM calculator (supporting the method family 'QMMM')
   */
  explicit QmmmGeometryOptimizer(Core::Calculator& calculator) : calculator_(calculator) {
    if (!calculator_.supportsMethodFamily("QMMM")) {
      throw std::runtime_error(
          "The calculator chosen for the QM/MM geometry optimization does not support the QM/MM method.");
    }
    auto silentLog = Core::Log::silent();
    calculator_.setLog(silentLog);
    fullOptimizer = std::make_unique<GeometryOptimizer<OptimizerType>>(calculator_);
    mmOptimizer = std::make_unique<GeometryOptimizer<OptimizerType>>(calculator_);
  };

  /**
   * @brief See GeometryOptimizerBase::optimize().
   *
   * @param atoms The molecular structure to be optimized.
   * @param log The logger.
   * @return int  The final number of full optimization cycles (for the full system) carried out.
   */
  int optimize(AtomCollection& atoms, Core::Log& log) {
    // Set settings for full optimizer
    Settings s1 = fullOptimizer->getSettings();
    s1.modifyString(SettingsNames::Optimizations::GeometryOptimizer::coordinateSystem,
                    CoordinateSystemInterpreter::getStringFromCoordinateSystem(this->coordinateSystem));
    s1.modifyInt(SettingsNames::Optimizations::Convergence::maxIter, this->maxFullOptMicrocycles);
    fullOptimizer->setSettings(s1);

    // Set settings for MM optimizer
    Settings s2 = mmOptimizer->getSettings();
    // Coordinate system must be Cartesian because of constraints
    s2.modifyString(SettingsNames::Optimizations::GeometryOptimizer::coordinateSystem,
                    CoordinateSystemInterpreter::getStringFromCoordinateSystem(CoordinateSystem::Cartesian));
    s2.modifyInt(SettingsNames::Optimizations::Convergence::maxIter, this->maxEnvOptMicrocycles);
    mmOptimizer->setSettings(s2);

    // Add atoms close to boundary to constrained atoms for optimization of the environment
    addBoundaryAtomsToConstrainedAtoms(atoms);

    // Perform optimization:
    bool converged = false;
    int cycles = 0;
    int totalCycles = 0;

    while (!converged) {
      cycles++;

      // MM-only opt of the environment
      if (cycles <= envOptSwitchOff) {
        if (envOptAtStart || cycles != 1) {
          calculator_.settings().modifyBool("ignore_qm", true);
          auto envMicrocycles = mmOptimizer->optimize(atoms, log);
          log.output << "Number of (MM-only) microcycles: " << envMicrocycles << Core::Log::endl;
          calculator_.settings().modifyBool("ignore_qm", false);
        }
      }

      // Full opt
      auto fullOptMicrocycles = fullOptimizer->optimize(atoms, log);
      if (fullOptMicrocycles < maxFullOptMicrocycles) {
        converged = true;
      }
      log.output << "Full optimization cycles in this macroiteration: " << fullOptMicrocycles << Core::Log::nl;
      totalCycles = (cycles - 1) * maxFullOptMicrocycles + fullOptMicrocycles;
      log.output << "Total full optimization cycles: " << totalCycles << Core::Log::endl;
      if (cycles == maxMacrocycles) {
        converged = true;
      }
    }
    return totalCycles;
  };

  /**
   * @brief Function to apply the given settings to underlying classes.
   * @param settings The new settings.
   */
  void setSettings(const Settings& settings) {
    fullOptimizer->check.applySettings(settings);
    mmOptimizer->check.applySettings(settings);
    fullOptimizer->optimizer.applySettings(settings);
    mmOptimizer->check.applySettings(settings);
    this->coordinateSystem = CoordinateSystemInterpreter::getCoordinateSystemFromString(
        settings.getString(SettingsNames::Optimizations::GeometryOptimizer::coordinateSystem));
    this->maxMacrocycles = settings.getInt(qmmmOptMaxMacroiterationsKey);
    this->maxFullOptMicrocycles = settings.getInt(qmmmOptMaxFullMicroiterationsKey);
    this->maxEnvOptMicrocycles = settings.getInt(qmmmOptMaxEnvMicroiterationsKey);
    this->boundaryDistanceThreshold =
        Utils::Constants::bohr_per_angstrom * settings.getDouble(qmmmOptBoundaryDistanceThresholdKey);
    this->envOptSwitchOff = settings.getInt(qmmmOptEnvOptSwitchOffKey);
    this->envOptAtStart = settings.getBool(qmmmOptEnvOptAtStartKey);
  };
  /**
   * @brief Get the public settings as a Utils::Settings object.
   * @return Settings A settings object with the current settings.
   */
  Settings getSettings() const {
    return QmmmGeometryOptimizerSettings<OptimizerType>(*this, *fullOptimizer, *mmOptimizer);
  };
  /**
   * @brief Get the settings of the calculator used for the energy calculations during the optimization.
   * @return std::shared_ptr<Settings> The settings of the calculator.
   */
  std::shared_ptr<Settings> getCalculatorSettings() const {
    return std::make_shared<Settings>(calculator_.settings());
  };
  /**
   * @brief Add an observer function that will be triggered in each iteration.
   *
   * @param function A function to be executed in every loop of the optimization.
   *                 The function will have access to the current cycle count,
   *                 the current value and to a const reference of the current
   *                 parameters.
   */
  void addObserver(std::function<void(const int&, const double&, const Eigen::VectorXd&)> function) {
    fullOptimizer->addObserver(function);
  }
  /**
   * @brief Clear all existing observer functions.
   *
   * For optimization problems with very fast evaluations of the underlying function
   * the removal of all observers can increase performance as the observers are given
   * as std::functions and can not be added via templates.
   */
  void clearObservers() {
    fullOptimizer->clearObservers();
    mmOptimizer->clearObservers();
  }
  /**
   * @brief Clears constrained atoms. Should be called between multiple optimize calls
   *
   */
  void clearConstrainedAtoms() {
    this->fullOptimizer->fixedAtoms = {};
    this->mmOptimizer->fixedAtoms = {};
  }
  /**
   * @brief The underlying geometry optimizer that optimizes the full system, public in order to change it's settings.
   */
  std::unique_ptr<GeometryOptimizer<OptimizerType>> fullOptimizer;
  /**
   * @brief The underlying geometry optimizer that optimizes the environment only,
   *        public in order to change it's settings.
   */
  std::unique_ptr<GeometryOptimizer<OptimizerType>> mmOptimizer;
  /// @brief The maximum number of macrocycles allowed (One cycle: Full opt. -> MM-only opt.).
  int maxMacrocycles = 30;
  /// @brief The maximum number of full optimization microcycles allowed per macrocycle.
  int maxFullOptMicrocycles = 15;
  /// @brief The maximum number of MM-only optimization microcycles allowed per macrocycle.
  int maxEnvOptMicrocycles = 1000;
  /**
   * @brief The distance threshold (in bohrs) determining which atoms of the environment
   *        are also frozen during the MM-only micro-opt.
   */
  double boundaryDistanceThreshold = 4.0 * Utils::Constants::bohr_per_angstrom;
  /// @brief The number of macrocycles after which the MM-only optimization is turned off.
  int envOptSwitchOff = 12;
  /// @brief Perform the MM-only optimization of the environment at the very beginning.
  bool envOptAtStart = true;
  /**
   * @brief Set the coordinate system in which the optimization shall be performed
   *
   * The optimization can be carried out in the internal coordinate space or with removed translations and rotations
   * possibly accelerating convergence.
   */
  CoordinateSystem coordinateSystem = CoordinateSystem::Internal;

 private:
  Core::Calculator& calculator_;
  // Function that adds the MM atoms close to the boundary to the constrained atoms during the MM-only opt.
  void addBoundaryAtomsToConstrainedAtoms(const Utils::AtomCollection& atoms) {
    auto mmOptSettings = mmOptimizer->getSettings();
    auto constrainedAtoms = calculator_.settings().getIntList("qm_atoms");

    // Add atoms close to boundary to the constrained atoms in MM optimization:
    std::vector<int> qmAtoms = mmOptimizer->fixedAtoms;
    // Loop over all atoms
    for (int i = 0; i < atoms.size(); ++i) {
      // Check whether it is not a QM atom
      if (std::find(qmAtoms.begin(), qmAtoms.end(), i) == qmAtoms.end()) {
        double minimalDistanceSquared = std::numeric_limits<double>::max();
        for (const auto& qmAtomIndex : qmAtoms) {
          double currentSquaredDistance = (atoms.getPosition(qmAtomIndex) - atoms.getPosition(i)).squaredNorm();
          if (currentSquaredDistance < minimalDistanceSquared) {
            minimalDistanceSquared = currentSquaredDistance;
          }
        }
        if (minimalDistanceSquared < std::pow(boundaryDistanceThreshold, 2)) {
          constrainedAtoms.push_back(i);
        }
      }
    }

    mmOptSettings.modifyIntList(SettingsNames::Optimizations::GeometryOptimizer::fixedAtoms, constrainedAtoms);
    mmOptimizer->setSettings(mmOptSettings);
  };
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_QMMMGEOMETRYOPTIMIZER_H
