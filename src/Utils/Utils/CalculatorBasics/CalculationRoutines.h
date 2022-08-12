/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_CALCULATION_ROUTINES_H_
#define UTILS_CALCULATION_ROUTINES_H_

#include <Core/Exceptions.h>
#include <Core/Interfaces/Calculator.h>
#include <Core/Log.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/Scf/LcaoUtils/SpinMode.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <boost/exception/diagnostic_information.hpp>
#include <string>

namespace Scine {
namespace Utils {
namespace CalculationRoutines {

/**
 * @brief Call calculator's calculate function with try catch and throw specific error for failure
 *
 * @param calculator A Scine calculator
 * @param log A Scine logger
 * @param errorMessage The error message of the exception
 * @throw Core::UnsuccessfulCalculationException if calculation fails
 * @return Results A Scine Results object
 */
inline Results calculateWithCatch(Core::Calculator& calculator, Core::Log& log, const std::string& errorMessage) {
  Utils::Results results;
  try {
    results = calculator.calculate("");
    if (!results.get<Property::SuccessfulCalculation>()) {
      throw Core::UnsuccessfulCalculationException("Calculator signalled unsuccessful calculation.");
    }
  }
  catch (...) {
    log.error << errorMessage << Core::Log::nl;
    throw;
  }
  return results;
}

/**
 * @brief Give calculator a new logger based on given booleans
 */
inline void setLog(Core::Calculator& calculator, bool error, bool warning, bool out) {
  auto log = Core::Log::silent();
  if (out) {
    log.output.add("cout", Core::Log::coutSink());
  }
  if (warning) {
    log.warning.add("cerr", Core::Log::cerrSink());
  }
  if (error) {
    log.error.add("cerr", Core::Log::cerrSink());
  }
  calculator.setLog(log);
}

/**
 * @brief Check all valid spin multiplicities within given range and change calculator to lowest energy spin multiplicity
 *
 * @param calculator A Scine calculator
 * @param log A Scine logger
 * @param maxSpinDifference The range of the multiplicity checks, e.g. 1 will check m=1 and m=5 for a m=3 system
 */
inline std::unique_ptr<Core::Calculator> spinPropensity(Core::Calculator& calculator, Core::Log& log, int maxSpinDifference) {
  // references
  auto referenceResults = calculator.results();
  auto refSettings = calculator.settings();
  int refSpinMultiplicity = refSettings.getInt(Utils::SettingsNames::spinMultiplicity);
  if (!referenceResults.has<Property::Energy>()) {
    calculator.setRequiredProperties(Property::Energy);
    referenceResults = calculator.calculate();
  }
  double refEnergy = referenceResults.get<Property::Energy>();

  // clone calculator with silent log to clone all calculators from this one
  auto bestClone = calculator.clone();
  setLog(*bestClone, false, false, false);
  double bestCloneEnergy = std::numeric_limits<double>::max();

  // set to unrestricted if it was restricted and singlet before
  if (refSpinMultiplicity == 1 && refSettings.getString(Utils::SettingsNames::spinMode) ==
                                      SpinModeInterpreter::getStringFromSpinMode(SpinMode::Restricted)) {
    bestClone->settings().modifyString(Utils::SettingsNames::spinMode,
                                       SpinModeInterpreter::getStringFromSpinMode(SpinMode::Unrestricted));
  }

  // determine lowest possible multiplicity above 0 and within range
  int lowestMultiplicity = std::numeric_limits<int>::max();
  for (int i = refSpinMultiplicity; i > 0 && i >= refSpinMultiplicity - (2 * maxSpinDifference); i = i - 2) {
    lowestMultiplicity = i;
  }

  // cycle all multiplicities within range
  for (int i = lowestMultiplicity; i <= refSpinMultiplicity + (2 * maxSpinDifference); i = i + 2) {
    if (i == refSpinMultiplicity) {
      continue;
    }
    auto clone = bestClone->clone();
    clone->settings().modifyInt(SettingsNames::spinMultiplicity, i);
    Results result;
    try {
      result = calculateWithCatch(*clone, log, "");
    }
    catch (...) {
      continue; // calculation failed, skip to next multiplicity
    }
    double propensityEnergy = result.get<Property::Energy>();
    if (propensityEnergy < bestCloneEnergy) {
      bestClone = clone->clone();
      bestCloneEnergy = propensityEnergy;
    }
  }
  if (bestCloneEnergy < refEnergy) {
    log.warning << "PropensityWarning: Detected spin multiplicity propensity; changing multiplicity to "
                << std::to_string(bestClone->settings().getInt(SettingsNames::spinMultiplicity)) << Core::Log::nl;
    return bestClone;
  }
  // no propensity
  return calculator.clone();
}

/**
 * @brief Split input string into method and dispersion string
 *
 * @param input a string specifying both with '-' as delimiter
 * @throw std::logic_error if more than one '-' present
 * @return std::pair of method and dispersion correction as strings
 */
inline std::pair<std::string, std::string> splitIntoMethodAndDispersion(const std::string& input) {
  if (input.empty()) {
    return std::make_pair("", "");
  }
  // check for exceptions
  // if input contains one of these as a substring
  // we return the input as method and empty dispersion
  std::vector<std::string> exceptions = {
      "PNO-CC",
      "HF-3C",
      "PBEH-3C",
      "B97-3C",
  };
  std::string inputCopy(input.size(), 0);
  std::transform(input.begin(), input.end(), inputCopy.begin(), [](unsigned char c) { return std::toupper(c); });
  if (std::any_of(exceptions.begin(), exceptions.end(),
                  [&, inputCopy](const std::string& e) { return inputCopy.find(e) != std::string::npos; }))
    return std::make_pair(input, "");

  std::string segment;
  std::vector<std::string> segments;
  std::stringstream ss(input);
  while (std::getline(ss, segment, '-')) {
    segments.push_back(segment);
  }
  if (segments.size() > 2) {
    throw std::logic_error("The provided method '" + input +
                           "' includes multiple '-' delimiter. Could not split into method and dispersion correction");
  }
  // we do not allow space in method specification to avoid unwanted calculation behavior
  if (segments[0].find(' ') != std::string::npos) {
    throw std::logic_error("The provided method '" + input + "' includes an empty space. This is currently not allowed.");
  }
  std::string dispersion = (segments.size() == 1) ? "" : segments[1];
  return std::make_pair(segments[0], dispersion);
}

inline void checkValidityOfChargeAndMultiplicity(int molecularCharge, int spinMultiplicity, const AtomCollection& atoms) {
  int numberUnpairedElectrons = spinMultiplicity - 1;

  int numElectrons = 0;
  for (const auto& atom : atoms) {
    numElectrons += Utils::ElementInfo::Z(atom.getElementType());
  }
  // Subtract the charge from the total number of electrons
  numElectrons -= molecularCharge;

  bool numUnpairedElecIsEven = numberUnpairedElectrons % 2 == 0;
  bool numElectronsIsEven = numElectrons % 2 == 0;
  if ((numUnpairedElecIsEven && !numElectronsIsEven) || (!numUnpairedElecIsEven && numElectronsIsEven)) {
    throw std::logic_error("Invalid charge/multiplicity pair for the given system!");
  }
}

} // namespace CalculationRoutines
} // namespace Utils
} // namespace Scine

#endif // UTILS_CALCULATION_ROUTINES_H_
