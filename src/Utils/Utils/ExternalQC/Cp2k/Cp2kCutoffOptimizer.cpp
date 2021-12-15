/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/ExternalQC/Cp2k/Cp2kCutoffOptimizer.h"
#include "Utils/CalculatorBasics/CalculationRoutines.h"
#include "Utils/ExternalQC/Cp2k/Cp2kCalculator.h"
#include "Utils/ExternalQC/Cp2k/Cp2kCalculatorSettings.h"
#include <iomanip>

namespace Scine {
namespace Utils {
namespace ExternalQC {

void Cp2kCutoffOptimizer::determineOptimalGridCutoffs(double energyEps, double distributionEpsFactor,
                                                      double startCutoff, double startRelCutoff) {
  if (distributionEpsFactor >= 1.0) {
    throw std::logic_error("The distribution threshold factor must be less than 1. Please read the documentation.");
  }
  if (!_calculator.getStructure()) {
    throw std::runtime_error("The given calculator does not hold a structure.");
  }
  if (_calculator.name() != "CP2K") {
    throw std::logic_error("The Cp2kCutoffOptimizer received a calculator that is not a Cp2kCalculator.");
  }
  _energyEps = energyEps;
  _distributionEpsFactor = distributionEpsFactor;
  // save existing settings to restore later
  auto prevSettings = Utils::Settings(_calculator.settings());
  // modify necessary settings for fast grid evaluation calculations
  _calculator.settings().modifyInt(Utils::SettingsNames::maxScfIterations, 1);
  _calculator.settings().modifyInt(SettingsNames::outerScf, 0);
  _calculator.settings().modifyBool(SettingsNames::allowUnconvergedScf, true);
  _calculator.setRequiredProperties(PropertyList(Property::Energy | Property::GridOccupation));
  if (_calculator.settings().getString(SettingsNames::scfGuess) == "restart") {
    _calculator.settings().modifyString(SettingsNames::scfGuess, "atomic");
  }
  std::pair<double, double> optimizedCutoffs = std::make_pair(startCutoff, startRelCutoff);
  try {
    // optimize cutoffs
    for (int i = 0; i < 3; ++i) {
      // repeat few times to be sure about convergence
      // this does not increase computation time heavily due to dataContainer
      auto cutoff = convergeCutoff(optimizedCutoffs.first, optimizedCutoffs.second);
      auto relCutoff = convergeCutoff(optimizedCutoffs.second, cutoff, true);
      optimizedCutoffs = convergeDistribution(cutoff, relCutoff);
    }
  }
  catch (...) {
    // ensure that original settings are restored
    _calculator.settings() = prevSettings;
    _calculator.getLog().error << "Failed to determine optimal grid cutoffs" << Core::Log::nl;
    throw;
  }
  // restore original settings
  _calculator.settings() = prevSettings;
  // set new cutoffs
  _calculator.settings().modifyDouble(SettingsNames::planeWaveCutoff, optimizedCutoffs.first);
  _calculator.settings().modifyDouble(SettingsNames::relMultiGridCutoff, optimizedCutoffs.second);
}

double Cp2kCutoffOptimizer::convergeCutoff(double startCutoff, double fixedCutoff, bool convergeRelCutoff) {
  double stepSize = (convergeRelCutoff) ? _relCutoffStepSize : _cutoffStepSize;
  double hardLimit = (convergeRelCutoff) ? _relCutoffHardLimit : _cutoffHardLimit;
  double maxCutoff = startCutoff; // upper start, may be increased in while loop
  if (convergeRelCutoff) {
    _calculator.settings().modifyDouble(SettingsNames::planeWaveCutoff, fixedCutoff);
  }
  else {
    _calculator.settings().modifyDouble(SettingsNames::relMultiGridCutoff, fixedCutoff);
  }
  while (true) {
    avoidInfiniteLoop(maxCutoff, hardLimit, fixedCutoff, convergeRelCutoff);
    auto refData = (convergeRelCutoff) ? getGridData(fixedCutoff, maxCutoff) : getGridData(maxCutoff, fixedCutoff);
    double refEnergy = refData.getEnergy();
    std::unique_ptr<double> oldCutoff;
    // cycle down
    for (double cutoff = maxCutoff - stepSize; cutoff > 0.0; cutoff -= stepSize) {
      auto data = (convergeRelCutoff) ? getGridData(fixedCutoff, cutoff) : getGridData(cutoff, fixedCutoff);
      double energy = data.getEnergy();
      if (std::fabs(energy - refEnergy) > _energyEps) {
        if (!oldCutoff) {        // this means the highest one was already outside of eps
          maxCutoff += stepSize; // need to calculate higher cutoff
          break;
        }
        else {
          // this means that we hit one that is outside, but the previous one is the lowest that is
          // still within the limit
          return *oldCutoff;
        }
      }
      // we are within eps -> can go lower
      oldCutoff = std::make_unique<double>(cutoff);
    }
    // if all values were within the limit, we return the lowest one
    if (oldCutoff && *oldCutoff < stepSize) {
      return *oldCutoff;
    }
  }
}

std::pair<double, double> Cp2kCutoffOptimizer::convergeDistribution(double cutoff, double relCutoff) {
  int nGrids = _dataContainer.getData(cutoff, relCutoff).getGridCounts().size();
  double idealDistribution = 1.0 / nGrids;
  double eps = _distributionEpsFactor * idealDistribution;
  while (true) {
    if (cutoff > _cutoffHardLimit || relCutoff > _relCutoffHardLimit) {
      throw std::runtime_error("The distribution could not be converged, stopped with cutoffs\n" +
                               std::to_string(cutoff) + " " + std::to_string(relCutoff));
    }
    auto data = getGridData(cutoff, relCutoff);
    auto candidateGridCounts = data.getGridCounts();
    int nGaussians = std::accumulate(candidateGridCounts.begin(), candidateGridCounts.end(),
                                     decltype(candidateGridCounts)::value_type(0));
    auto minPtr = min_element(candidateGridCounts.begin(), candidateGridCounts.end());
    if (*minPtr / static_cast<double>(nGaussians) < eps) { // not converged, check which to change
      if (minPtr - candidateGridCounts.begin() > nGrids / 2.0) {
        cutoff += _cutoffStepSize;
      }
      else {
        relCutoff += _relCutoffStepSize;
      }
    }
    else {
      return std::make_pair(cutoff, relCutoff);
    }
  }
}

Cp2kCutoffData Cp2kCutoffOptimizer::getGridData(double cutoff, double relCutoff) {
  if (_dataContainer.has(cutoff, relCutoff)) {
    return _dataContainer.getData(cutoff, relCutoff);
  }
  _calculator.settings().modifyDouble(SettingsNames::planeWaveCutoff, cutoff);
  _calculator.settings().modifyDouble(SettingsNames::relMultiGridCutoff, relCutoff);
  Results results = CalculationRoutines::calculateWithCatch(_calculator, _calculator.getLog(),
                                                            "Aborting cutoff optimization due to failed calculation "
                                                            "with cutoffs " +
                                                                std::to_string(cutoff) + " and " + std::to_string(relCutoff));
  auto gridCounts = results.get<Property::GridOccupation>();
  double energy = results.get<Property::Energy>();
  auto newData = ExternalQC::Cp2kCutoffData(cutoff, relCutoff, energy, gridCounts);
  _dataContainer.add(newData);
  return newData;
}

void Cp2kCutoffOptimizer::avoidInfiniteLoop(double maxCutoff, double hardLimit, double fixedCutoff, bool convergeRelCutoff) const {
  if (maxCutoff > hardLimit) {
    if (convergeRelCutoff) {
      throw std::runtime_error("Reached the maximum relCutoff of " + std::to_string(hardLimit) +
                               " (with the cutoff "
                               "previously converged to " +
                               std::to_string(fixedCutoff) +
                               ") without finding convergence.\n"
                               "Change your start cutoffs or the energy threshold.");
    }
    throw std::runtime_error("Reached the maximum cutoff of " + std::to_string(hardLimit) +
                             " with the relative cutoff "
                             "fixed to " +
                             std::to_string(fixedCutoff) +
                             " without finding convergence.\n"
                             "Change your start cutoffs or the energy threshold.");
  }
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
