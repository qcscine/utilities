/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILSOS_DAVIDSONSETTINGS_H
#define UTILSOS_DAVIDSONSETTINGS_H

#include "SubspaceCollapser.h"
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>

namespace Scine {
namespace Utils {

class InvalidDiagonalizerInputException : public std::exception {
 public:
  explicit InvalidDiagonalizerInputException(std::string error) : error_("Input error: " + error) {
  }
  const char* what() const noexcept final {
    return error_.c_str();
  }

 private:
  std::string error_;
};

static constexpr const char* numberOfRootsOption = "number_of_roots";
static constexpr const char* initialGuessDimensionOption = "initial_guess_dimension";
static constexpr const char* seedOption = "seed";
static constexpr const char* residualNormToleranceOption = "residual_norm_tolerance";
static constexpr const char* correctionToleranceOption = "correction_tolerance";
static constexpr const char* subspaceCollapseDimensionOption = "collapse_dimension";

static constexpr const char* gepAlgorithmForBalancedMethodOption = "gep_algo";

/**
 * @brief Settings class for the iterative diagonalizers.
 */
class DiagonalizerSettings : public Utils::Settings {
 public:
  DiagonalizerSettings(int eigenvaluesToCompute, int totalDimension) : Utils::Settings("Diagonalizer settings") {
    UniversalSettings::IntDescriptor numberOfRoots("Number of roots to diagonalize");
    numberOfRoots.setMinimum(1);
    numberOfRoots.setDefaultValue(eigenvaluesToCompute);
    _fields.push_back(numberOfRootsOption, std::move(numberOfRoots));

    UniversalSettings::IntDescriptor initialGuessDimension("Number of initial guess vectors");
    initialGuessDimension.setMinimum(eigenvaluesToCompute);
    initialGuessDimension.setDefaultValue(eigenvaluesToCompute);
    _fields.push_back(initialGuessDimensionOption, std::move(initialGuessDimension));

    UniversalSettings::IntDescriptor maxDavidsonIterations("Number of maximal iterations");
    maxDavidsonIterations.setMinimum(1);
    maxDavidsonIterations.setMaximum(totalDimension);
    maxDavidsonIterations.setDefaultValue(totalDimension);
    _fields.push_back(SettingsNames::maxDavidsonIterations, std::move(maxDavidsonIterations));

    UniversalSettings::IntDescriptor seed("Seed for the random guess initializer");
    seed.setDefaultValue(42);
    _fields.push_back(seedOption, std::move(seed));

    UniversalSettings::DoubleDescriptor eigenvalueTolerance("Convergence threshold for the eigenvalues");
    eigenvalueTolerance.setMinimum(std::numeric_limits<double>::min());
    eigenvalueTolerance.setDefaultValue(1e-5);
    _fields.push_back(residualNormToleranceOption, std::move(eigenvalueTolerance));

    resetToDefaults();
    check(totalDimension);
  }

  int calculateInitialSubspaceCollapserIterations(int eigenvaluesToCompute, int maxDimension) const {
    return SubspaceCollapser::calculateSubspaceCollapserIterations(eigenvaluesToCompute, 0, maxDimension);
  }

  void check(int totalDimension) {
    int eigenvaluesToCompute = getInt(numberOfRootsOption);
    int subspaceDimension = getInt(initialGuessDimensionOption);

    if (eigenvaluesToCompute < 0) {
      throw InvalidDiagonalizerInputException(
          "Unintended behaviour: calculate negative amount of eigenvalues in diagonalizer.");
    }

    if (eigenvaluesToCompute > totalDimension) {
      throw InvalidDiagonalizerInputException("Number of eigenvalues sought is higher than the total dimension.");
    }
    if (subspaceDimension < eigenvaluesToCompute || subspaceDimension > totalDimension) {
      throw InvalidDiagonalizerInputException(
          "Subspace dimension initially smaller than the number of eigenvalues to compute"
          "or subspace dimension bigger than the total dimension.");
    }
    modifyInt(initialGuessDimensionOption, subspaceDimension);
  }
};

/**
 * @brief Class for the settings of digonalizers based on Krylov methods, i.e. Davidson.
 * @class KrylovSettings
 */
class KrylovSettings : public DiagonalizerSettings {
 public:
  KrylovSettings(int eigenvaluesToCompute, int totalDimension)
    : DiagonalizerSettings(eigenvaluesToCompute, totalDimension) {
    UniversalSettings::DoubleDescriptor correctionTolerance("Acceptance threshold for correction vectors");
    correctionTolerance.setMinimum(std::numeric_limits<double>::min());
    correctionTolerance.setDefaultValue(5e-4);
    _fields.push_back(correctionToleranceOption, std::move(correctionTolerance));

    UniversalSettings::IntDescriptor subspaceCollapseDimension(
        "Maximal dimension of the subspace after which to collapse");
    subspaceCollapseDimension.setMinimum(2 * eigenvaluesToCompute);
    subspaceCollapseDimension.setDefaultValue(2 * eigenvaluesToCompute);
    _fields.push_back(subspaceCollapseDimensionOption, std::move(subspaceCollapseDimension));

    UniversalSettings::OptionListDescriptor gepAlgorithm("Algorithm to compute the stable Generalized"
                                                         "Eigenvalue Problem Ax=lBx when B is almost singular.");
    gepAlgorithm.addOption("standard");
    gepAlgorithm.addOption("cholesky");
    gepAlgorithm.addOption("simultaneous_diag");
    gepAlgorithm.setDefaultOption("simultaneous_diag");
    _fields.push_back(gepAlgorithmForBalancedMethodOption, std::move(gepAlgorithm));

    resetToDefaults();
    check(totalDimension);
    modifyInt(subspaceCollapseDimensionOption,
              calculateInitialSubspaceCollapserIterations(getInt(numberOfRootsOption), totalDimension));
  }
};

} // namespace Utils
} // namespace Scine

#endif // UTILSOS_DAVIDSONSETTINGS_H
