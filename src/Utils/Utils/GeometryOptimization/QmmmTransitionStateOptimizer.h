/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_QMMMTRANSITIONSTATEOPTIMIZER_H
#define UTILS_QMMMTRANSITIONSTATEOPTIMIZER_H

#include "GeometryOptimizer.h"
#include <Core/Interfaces/EmbeddingCalculator.h>
#include <Core/Log.h>
#include <Utils/Geometry.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Typenames.h>
#include <Utils/UniversalSettings/OptimizationSettingsNames.h>
#include <Utils/UniversalSettings/SettingsNames.h>

namespace Scine {
namespace Utils {
/* Forward declaration */
template<class OptimizerType, class MmOptimizerType>
class QmmmTransitionStateOptimizer;

/**
 * @class QmmmTransitionStateOptimizerSettings QmmmTransitionStateOptimizer.h
 * @brief Settings for a QmmmTransitionStateOptimizer
 *
 * Uses template arguments in order to automatically include the
 * settings of underlying objects into the given settings.
 *
 * @tparam OptimizerType The underlying Optimizer class.
 */
template<class OptimizerType, class MmOptimizerType>
class QmmmTransitionStateOptimizerSettings : public Settings {
 public:
  /**
   * @brief Constructor.
   *
   * @param qmmmTSOptimizer The QM/MM geometry optimizer.
   * @param qmOptimizer The geometry optimizer for the full system.
   * @param mmOptimizer The geometry optimizer for the environment.
   */
  QmmmTransitionStateOptimizerSettings(const QmmmTransitionStateOptimizer<OptimizerType, MmOptimizerType>& qmmmTSOptimizer,
                                       const GeometryOptimizer<OptimizerType>& qmOptimizer,
                                       const GeometryOptimizer<MmOptimizerType>& mmOptimizer)
    : Settings("QmmmTransitionStateOptimizerSettings") {
    /*
     * The full system optimizer and the environment-only optimizers always have the same convergence criteria.
     *
     * The setting "max_iter" is overwritten by the special settings regarding the maximum number of iterations
     * that are present in the QM/MM optimizer settings. It is therefore meaningless to modify that setting
     * in the convergence check.
     */
    qmOptimizer.optimizer.addSettingsDescriptors(this->_fields);
    qmOptimizer.check.addSettingsDescriptors(this->_fields);
    mmOptimizer.optimizer.addSettingsDescriptors(this->_fields);

    UniversalSettings::OptionListDescriptor geooptCoordinateSystem("Set the coordinate system.");
    geooptCoordinateSystem.addOption("internal");
    geooptCoordinateSystem.addOption("cartesianWithoutRotTrans");
    geooptCoordinateSystem.addOption("cartesian");
    geooptCoordinateSystem.setDefaultOption(
        CoordinateSystemInterpreter::getStringFromCoordinateSystem(qmmmTSOptimizer.coordinateSystem));
    this->_fields.push_back(SettingsNames::Optimizations::GeometryOptimizer::coordinateSystem,
                            std::move(geooptCoordinateSystem));

    UniversalSettings::IntDescriptor qmmmOptMaxMacroiterations("The maximum number of macrocycles allowed.");
    qmmmOptMaxMacroiterations.setDefaultValue(qmmmTSOptimizer.maxMacrocycles);
    this->_fields.push_back(QmmmTransitionStateOptimizer<OptimizerType, MmOptimizerType>::qmmmOptMaxMacroiterationsKey,
                            std::move(qmmmOptMaxMacroiterations));

    UniversalSettings::IntDescriptor qmmmOptMaxEnvMicroiterations(
        "The maximum number of MM-only optimization microcycles allowed per macrocycle.");
    qmmmOptMaxEnvMicroiterations.setDefaultValue(qmmmTSOptimizer.maxEnvOptMicrocycles);
    this->_fields.push_back(QmmmTransitionStateOptimizer<OptimizerType, MmOptimizerType>::qmmmOptMaxEnvMicroiterationsKey,
                            std::move(qmmmOptMaxEnvMicroiterations));

    UniversalSettings::IntDescriptor qmmmOptMaxQmMicroiterations(
        "The maximum number of QM optimization microcycles allowed per macrocycle.");
    qmmmOptMaxQmMicroiterations.setDefaultValue(qmmmTSOptimizer.maxQmOptMicrocycles);
    this->_fields.push_back(QmmmTransitionStateOptimizer<OptimizerType, MmOptimizerType>::qmmmOptMaxQmMicroiterationsKey,
                            std::move(qmmmOptMaxQmMicroiterations));

    // keep
    UniversalSettings::IntDescriptor qmmmOptEnvOptSwitchOff(
        "The number of macrocycles after which the MM-only optimization is switched off.");
    qmmmOptEnvOptSwitchOff.setDefaultValue(qmmmTSOptimizer.envOptSwitchOff);
    this->_fields.push_back(QmmmTransitionStateOptimizer<OptimizerType, MmOptimizerType>::qmmmOptEnvOptSwitchOffKey,
                            std::move(qmmmOptEnvOptSwitchOff));

    // keep
    UniversalSettings::BoolDescriptor qmmmOptEnvOptAtStart(
        "Perform the MM-only optimization of the environment at the very beginning.");
    qmmmOptEnvOptAtStart.setDefaultValue(qmmmTSOptimizer.envOptAtStart);
    this->_fields.push_back(QmmmTransitionStateOptimizer<OptimizerType, MmOptimizerType>::qmmmOptEnvOptAtStartKey,
                            std::move(qmmmOptEnvOptAtStart));
    this->resetToDefaults();
  }
};

/**
 * @class QmmmTransitionStateOptimizer QmmmTransitionStateOptimizer.h
 *
 * @brief Implementation of a QM/MM geometry optimizer that should be preferably applied for QM/MM geometry
 *        optimizations as it employs the microiteration-optimization scheme leading to a more efficient optimization
 *        with regards to the number of QM gradients calculations needed.
 *
 * Basic algorithm:
 * 1. optimize only the environment (MM region) until convergence or maximum number of
 *    iterations reached (setting: qmmm_opt_max_env_microiterations).
 * 2. optimize only the QM system with the QM method until convergence or maximum number of
 *    iterations reached (setting: qmmm_opt_max_qm_microiterations).
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
 * @tparam MmOptimizerType Expects any of the Optimizer classes. Note that some special optimizers
 *                       may not yet be supported or may need additional specialization.
 */
template<class OptimizerType, class MmOptimizerType>
class QmmmTransitionStateOptimizer : public GeometryOptimizerBase {
 public:
  static constexpr const char* qmmmOptMaxMacroiterationsKey = "qmmm_opt_max_macroiterations";
  static constexpr const char* qmmmOptMaxEnvMicroiterationsKey = "qmmm_opt_max_env_microiterations";
  static constexpr const char* qmmmOptMaxQmMicroiterationsKey = "qmmm_opt_max_qm_microiterations";
  static constexpr const char* qmmmOptEnvOptAtStartKey = "qmmm_opt_env_start";
  static constexpr const char* qmmmOptEnvOptSwitchOffKey = "qmmm_opt_env_switch_off";
  /**
   * @brief Constructor.
   * @param calculator The calculator to be used for the single point/gradient calculations.
   *                   The calculator has to be a QM/MM calculator (supporting the method family 'QMMM')
   */
  explicit QmmmTransitionStateOptimizer(std::shared_ptr<Core::Calculator>& qmmmCalculator)
    : qmmmCalculator_(qmmmCalculator) {
    if (!qmmmCalculator_->supportsMethodFamily("QMMM")) {
      throw std::runtime_error(
          "The calculator chosen for the QM/MM geometry optimization does not support the QM/MM method.");
    }
    // cast the calculator to an embedding calculator
    auto castedCalc = std::dynamic_pointer_cast<Scine::Core::EmbeddingCalculator>(qmmmCalculator);
    if (!castedCalc) {
      throw std::runtime_error("Please specify an embedding calculator.");
    }
    auto underlyingCalculators = castedCalc->getUnderlyingCalculators();

    auto silentLog = Core::Log::silent();
    qmmmCalculator_->setLog(silentLog);
    qmCalculator_ = underlyingCalculators[0]->clone();
    auto mmCalculator = underlyingCalculators[1];
    qmOptimizer = std::make_shared<GeometryOptimizer<OptimizerType>>(*qmCalculator_);
    mmOptimizer = std::make_shared<GeometryOptimizer<MmOptimizerType>>(*qmmmCalculator_);
  };

  /**
   * @brief See GeometryOptimizerBase::optimize().
   *
   * @param atoms The molecular structure to be optimized.
   * @param log The logger.
   * @return int  The final number of full optimization cycles (for the full system) carried out.
   */
  int optimize(AtomCollection& atoms, Core::Log& log) override {
    qmmmCalculator_->setStructure(atoms);
    auto castedCalc = std::dynamic_pointer_cast<Scine::Core::EmbeddingCalculator>(qmmmCalculator_);
    if (!castedCalc) {
      throw std::runtime_error("Please specify an embedding calculator.");
    }
    auto underlyingCalculators = castedCalc->getUnderlyingCalculators();
    qmCalculator_->setStructure(*underlyingCalculators[0]->getStructure());
    // Set settings for full optimizer
    Settings s1 = qmOptimizer->getSettings();
    s1.modifyString(SettingsNames::Optimizations::GeometryOptimizer::coordinateSystem,
                    CoordinateSystemInterpreter::getStringFromCoordinateSystem(this->coordinateSystem));
    s1.modifyInt(SettingsNames::Optimizations::Convergence::maxIter, this->maxQmOptMicrocycles);
    qmOptimizer->setSettings(s1);
    // Set settings for MM optimizer
    Settings s2 = mmOptimizer->getSettings();
    // Coordinate system must be Cartesian because of constraints
    s2.modifyString(SettingsNames::Optimizations::GeometryOptimizer::coordinateSystem,
                    CoordinateSystemInterpreter::getStringFromCoordinateSystem(CoordinateSystem::Cartesian));
    s2.modifyInt(SettingsNames::Optimizations::Convergence::maxIter, this->maxEnvOptMicrocycles);
    mmOptimizer->setSettings(s2);

    // Obtain shift caused by optimization of QM region in non-Cartesian coordinate system
    if (qmOptimizer->coordinateSystem == CoordinateSystem::Internal) {
      Utils::PositionCollection updatedPositions = getShiftedMmAndQmPositions(atoms.getPositions(), false);
      // Updated relevant objects with shifted positions
      atoms.setPositions(updatedPositions);
      qmmmCalculator_->modifyPositions(updatedPositions);
      qmCalculator_->modifyPositions(underlyingCalculators[0]->getStructure()->getPositions());
    }
    else if (qmOptimizer->coordinateSystem == CoordinateSystem::CartesianWithoutRotTrans) {
      Utils::PositionCollection updatedPositions = getShiftedMmAndQmPositions(atoms.getPositions(), true);
      // Updated relevant objects with shifted positions
      atoms.setPositions(updatedPositions);
      qmmmCalculator_->modifyPositions(updatedPositions);
      qmCalculator_->modifyPositions(underlyingCalculators[0]->getStructure()->getPositions());
    }
    else if (qmOptimizer->coordinateSystem == CoordinateSystem::Cartesian) {
      log.warning << "Warning: You chose a true Cartesian coordinate system for the QM optimization in a QM/MM "
                     "transition state optimization. This may cause constant translations during the optimizations "
                     "which can cause clashes between the optimized QM and frozen MM region"
                  << Core::Log::endl;
    }

    // Add QM atoms to constrained atoms for optimization of the environment
    addQmAtomsToConstrainedAtomsForMM();

    bool electrostaticEmbedding = qmmmCalculator_->settings().getBool(SettingsNames::electrostaticEmbedding);

    // Perform optimization:
    bool qmConverged = false;
    bool mmConverged = false;
    bool converged = false;
    int cycles = 0;
    int totalCycles = 0;
    double qmmmEnergy = 0.0;
    double mmEnergy = 0.0;
    auto qmAtomsList = qmmmCalculator_->settings().getIntList(SettingsNames::qmAtomsList);
    // handle observers, make sure ours is first
    auto defaultObservers = qmOptimizer->optimizer.getObservers();
    qmOptimizer->clearObservers();
    auto customObserver = [&](const int& /* cycle */, const double& /*energy*/, const Eigen::VectorXd& /*params*/) {
      updatePositionsInFullSystem(atoms, *qmCalculator_->getStructure(), qmAtomsList);
      qmmmCalculator_->modifyPositions(atoms.getPositions());
    };
    qmOptimizer->addObserver(customObserver);
    for (auto& observer : defaultObservers) {
      qmOptimizer->addObserver(observer);
    }
    // optimize until convergence
    while (!converged) {
      cycles++;
      // MM-only opt of the environment
      if (qmConverged || envOptAtStart || cycles != 1) {
        // optimize MM until convergence
        int totalMmCycles = 0;
        while (!mmConverged) {
          // if electrostatic embedding is applied, let the QM/MM calculator feel the MM charges
          qmmmCalculator_->settings().modifyBool(SettingsNames::ignoreQmOption, true);
          qmmmCalculator_->settings().modifyBool(SettingsNames::electrostaticEmbedding, false);
          qmmmCalculator_->setRequiredProperties(Property::Energy | Property::Gradients);
          qmmmCalculator_->setStructure(atoms);
          const auto preResults = qmmmCalculator_->calculate();
          const auto preEnergy = preResults.get<Property::Energy>();
          const GradientCollection preGrad = preResults.get<Property::Gradients>();
          auto check = mmOptimizer->getConvergenceCheck();
          const PositionCollection positions = qmmmCalculator_->getStructure()->getPositions();
          const Eigen::VectorXd param = Eigen::Map<const Eigen::VectorXd>(positions.data(), positions.size());
          check.setParametersAndValue(param, mmEnergy);

          bool haveToChangeMM =
              std::fabs(preEnergy - mmEnergy) >= check.deltaValue ||
              !check.checkConvergence(param, preEnergy, Eigen::Map<const Eigen::VectorXd>(preGrad.data(), preGrad.size()));
          int envMicroCycles = 0;
          if (haveToChangeMM) {
            envMicroCycles = mmOptimizer->optimize(atoms, log);
            totalMmCycles += envMicroCycles;
          }
          else {
            log.output << "# MM is already converged, skipping microcycles: " << totalMmCycles << Core::Log::endl;
          }
          // restore settings
          mmEnergy = preEnergy;
          qmmmCalculator_->settings().modifyBool(SettingsNames::ignoreQmOption, false);
          qmmmCalculator_->settings().modifyBool(SettingsNames::electrostaticEmbedding, electrostaticEmbedding);
          qmmmCalculator_->setRequiredProperties(Property::Energy | Property::Gradients); // required to receive
                                                                                          // correct properties
          qmmmCalculator_->results() = Results{}; // purge results to ensure we get QM
          if (haveToChangeMM) {
            // update positions
            qmmmCalculator_->modifyPositions(atoms.getPositions());
            // make sure that the QM calculator gets the new positions, but make them reasonable without Hessian first
            const bool optLinks = qmmmCalculator_->settings().getBool(SettingsNames::optimizeLinks);
            qmmmCalculator_->settings().modifyBool(SettingsNames::optimizeLinks, true);
            qmmmCalculator_->calculate();
            qmmmCalculator_->settings().modifyBool(SettingsNames::optimizeLinks, optLinks);
            qmCalculator_->modifyPositions(underlyingCalculators[0]->getStructure()->getPositions());
            log.output << "# MM-only microcycles: " << totalMmCycles << Core::Log::endl;
          }
          if (envMicroCycles < maxEnvOptMicrocycles) {
            mmConverged = true;
          }
        }
        // update the point charge positions for the QM calculator
        if (electrostaticEmbedding) {
          // make sure that the settings responsible for the point charges are set correctly
          qmCalculator_->settings() = underlyingCalculators[0]->settings();
          // make sure we get the right properties
          auto requiredProperties = qmCalculator_->getRequiredProperties();
          requiredProperties.addProperty(Property::PointChargesGradients);
          qmCalculator_->setRequiredProperties(requiredProperties);
        }
      }

      // Dimer specific settings to save senseless rotations
      if (cycles != 1 && qmOptimizer->getSettings().valueExists(SettingsNames::Optimizations::Dimer::skipFirstRotation)) {
        qmOptimizer->getSettings().modifyBool(SettingsNames::Optimizations::Dimer::skipFirstRotation, true);
      }
      auto qmStructure = qmCalculator_->getStructure();
      int qmOptMicrocycles = -1;
      // Catch exception in QM optimization
      try {
        qmOptMicrocycles = qmOptimizer->optimize(*qmStructure, log);
      }
      catch (const std::exception& error) {
        log.error << "QM-only optimization failed with:\n" << error.what() << Core::Log::endl;
        throw;
      }
      const double currentQmmmEnergy = (qmmmCalculator_->results().has<Property::Energy>())
                                           ? qmmmCalculator_->results().get<Property::Energy>()
                                           : qmmmCalculator_->calculate().get<Property::Energy>();
      // Wait for convergence of QM optimization
      if (qmOptMicrocycles < maxQmOptMicrocycles and qmOptMicrocycles != -1) {
        if (qmConverged) {
          if (std::abs(currentQmmmEnergy - qmmmEnergy) < getConvergenceCheck().deltaValue) {
            // basically convergence reached, breaking while loop for all systems
            break;
          }
        }
        qmmmEnergy = currentQmmmEnergy;
        qmConverged = true;
        mmConverged = false;
      }
      else {
        qmmmEnergy = currentQmmmEnergy;
        qmConverged = false;
        mmConverged = false;
      }

      log.output << "# QM-only optimization macrocycles: " << cycles << Core::Log::endl;
      if (qmOptMicrocycles != -1) {
        totalCycles += qmOptMicrocycles;
      }

      if (cycles == maxMacrocycles) {
        converged = true;
      }
    }
    log.output << "# QM-only optimization total cycles: " << totalCycles << Core::Log::endl;
    // reset settings of QM/MM calculator, because they have been changed internally many times
    return cycles;
  };

  /**
   * @brief Function to apply the given settings to underlying classes.
   * @param settings The new settings.
   */
  void setSettings(const Settings& settings) override {
    // apply settings to underlying optimizer
    qmOptimizer->check.applySettings(settings);
    mmOptimizer->check.applySettings(settings);
    mmOptimizer->optimizer.applySettings(settings);
    qmOptimizer->optimizer.applySettings(settings);
    this->coordinateSystem = CoordinateSystemInterpreter::getCoordinateSystemFromString(
        settings.getString(SettingsNames::Optimizations::GeometryOptimizer::coordinateSystem));
    // apply custom QMMM transition state optimization settings
    this->maxMacrocycles = settings.getInt(qmmmOptMaxMacroiterationsKey);
    this->maxEnvOptMicrocycles = settings.getInt(qmmmOptMaxEnvMicroiterationsKey);
    this->envOptAtStart = settings.getBool(qmmmOptEnvOptAtStartKey);
    this->maxQmOptMicrocycles = settings.getInt(qmmmOptMaxQmMicroiterationsKey);
    this->envOptSwitchOff = settings.getInt(qmmmOptEnvOptSwitchOffKey);
  };
  /**
   * @brief Get the public settings as a Utils::Settings object.
   * @return Settings A settings object with the current settings.
   */
  Settings getSettings() const override {
    return QmmmTransitionStateOptimizerSettings<OptimizerType, MmOptimizerType>(*this, *qmOptimizer, *mmOptimizer);
  };
  /**
   * @brief Get the settings of the calculator used for the energy calculations during the optimization.
   * @return std::shared_ptr<Settings> The settings of the calculator.
   */
  std::shared_ptr<Settings> getCalculatorSettings() const override {
    return std::make_shared<Settings>(qmmmCalculator_->settings());
  };
  /**
   * @brief Add an observer function that will be triggered in each iteration.
   *
   * @param function A function to be executed in every loop of the optimization.
   *                 The function will have access to the current cycle count,
   *                 the current value and to a const reference of the current
   *                 parameters.
   */
  void addObserver(std::function<void(const int&, const double&, const Eigen::VectorXd&)> function) override {
    qmOptimizer->addObserver(function);
  };
  /**
   * @brief Clear all existing observer functions.
   *
   * For optimization problems with very fast evaluations of the underlying function
   * the removal of all observers can increase performance as the observers are given
   * as std::functions and can not be added via templates.
   */
  void clearObservers() override {
    qmOptimizer->clearObservers();
    mmOptimizer->clearObservers();
  };
  /**
   * @brief Clears constrained atoms. Should be called between multiple optimize calls
   *
   */
  inline void clearConstrainedAtoms() {
    this->qmOptimizer->fixedAtoms = {};
  };

  inline void updatePositionsInFullSystem(Utils::AtomCollection& fullSystem, const Utils::AtomCollection& coreRegion,
                                          const std::vector<int> coreIndices) {
    int coreIndexCounter = 0;
    for (auto& id : coreIndices) {
      if (fullSystem.at(id).getElementType() != coreRegion.at(coreIndexCounter).getElementType()) {
        throw std::runtime_error("Cannot update position of QM region atoms in the full system!");
      }
      fullSystem.setPosition(id, coreRegion.getPosition(coreIndexCounter));
      coreIndexCounter++;
    }
  }

  const GradientBasedCheck& getConvergenceCheck() const override {
    return qmOptimizer->check;
  }
  /// @brief Virtual default destructor.
  virtual ~QmmmTransitionStateOptimizer() = default;

  /**
   * @brief The underlying geometry optimizer that optimizes the QM system, public in order to change it's settings.
   */
  std::shared_ptr<GeometryOptimizer<OptimizerType>> qmOptimizer;
  /**
   * @brief The underlying geometry optimizer that optimizes the environment only,
   *        public in order to change it's settings.
   */
  std::shared_ptr<GeometryOptimizer<MmOptimizerType>> mmOptimizer;
  /// @brief The maximum number of macrocycles allowed (One cycle: Full opt. -> MM-only opt.).
  int maxMacrocycles = 150;
  /// @brief The maximum number of full optimization microcycles allowed per macrocycle. The
  /// default is set to 14, because the mostly used optimizer (Bofill) has the default to calculate a Hessian every
  /// 5th step
  int maxQmOptMicrocycles = 14;
  /// @brief The maximum number of MM-only optimization microcycles allowed per macrocycle.
  int maxEnvOptMicrocycles = 1000;
  /// @brief The number of macrocycles after which the MM-only optimization is turned off.
  int envOptSwitchOff = 12;
  /// @brief Perform the MM-only optimization of the environment at the very beginning.
  bool envOptAtStart = true;
  /**
   * @brief Set the coordinate system in which the QM optimization shall be performed
   *
   * The optimization can be carried out in the internal coordinate space or with removed translations and rotations
   * possibly accelerating convergence.
   */
  CoordinateSystem coordinateSystem = CoordinateSystem::CartesianWithoutRotTrans;

 private:
  std::shared_ptr<Core::Calculator>& qmmmCalculator_;
  std::shared_ptr<Core::Calculator> qmCalculator_;
  std::vector<double> pcMatrix_;

  inline void addQmAtomsToConstrainedAtomsForMM() {
    auto mmOptSettings = mmOptimizer->getSettings();

    auto qmAtomsList = qmmmCalculator_->settings().getIntList(SettingsNames::qmAtomsList);

    mmOptSettings.modifyIntList(SettingsNames::Optimizations::GeometryOptimizer::fixedAtoms, qmAtomsList);
    mmOptimizer->setSettings(mmOptSettings);
  };

  inline Utils::PositionCollection getShiftedMmAndQmPositions(const Utils::PositionCollection& orgPositions, bool rotTransOnly) {
    // Mimik transformation of QM region to obtain shift vector
    Utils::AtomCollection tmpQmStructure = *qmCalculator_->getStructure();
    Utils::PositionCollection orgQmPositions = tmpQmStructure.getPositions();
    // Simulate transformation once to obtain shift vector and shift
    std::shared_ptr<InternalCoordinates> transformation = std::make_shared<InternalCoordinates>(tmpQmStructure, rotTransOnly);
    Eigen::VectorXd transformedPositions = transformation->coordinatesToInternal(orgQmPositions);
    Utils::PositionCollection secondTransShift = transformation->coordinatesToCartesian(transformedPositions);
    // Translation vector for COMs
    const auto masses = Utils::Geometry::Properties::getMasses(tmpQmStructure.getElements());
    const Eigen::RowVector3d& shift = Utils::Geometry::Properties::getCenterOfMass(secondTransShift, masses) -
                                      Utils::Geometry::Properties::getCenterOfMass(orgQmPositions, masses);
    // Updated all positions by change to come
    Utils::PositionCollection shiftedPositions = Utils::Geometry::Manipulations::translatePositions(orgPositions, shift);

    return shiftedPositions;
  }
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_QMMMTRANSITIONSTATEOPTIMIZER_H
