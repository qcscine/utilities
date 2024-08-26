/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_CALCULATION_ROUTINES_H_
#define UTILS_CALCULATION_ROUTINES_H_

#include <Core/Exceptions.h>
#include <Core/Interfaces/Calculator.h>
#include <Core/Interfaces/EmbeddingCalculator.h>
#include <Core/Log.h>
#include <Core/ModuleManager.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/Scf/LcaoUtils/SpinMode.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <boost/algorithm/string.hpp>
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

inline void inputPreparation(std::string& method_family, std::string& program) {
  if (!program.empty()) {
    std::transform(std::begin(program) + 1, std::end(program), std::begin(program) + 1,
                   [](unsigned char x) { return std::tolower(x); });

    program[0] = std::toupper(program[0]);
  }

  std::transform(std::begin(method_family), std::end(method_family), std::begin(method_family),
                 [](unsigned char x) { return std::toupper(x); });
}

inline std::vector<std::string> split(const std::string& s, char delimiter) {
  std::vector<std::string> elements;
  std::stringstream ss;
  ss.str(s);
  std::string item;

  while (std::getline(ss, item, delimiter)) {
    elements.push_back(item);
  }
  return elements;
}

/**
 * @brief Join a vector of strings with '/' as delimiter.
 * @param sVector    The strings.
 * @param joinFirstN Join only the first N elements in sVector.
 * @return The joined string.
 */
inline std::string join(const std::vector<std::string>& sVector, unsigned int joinFirstN) {
  std::string result;
  joinFirstN = (joinFirstN <= sVector.size()) ? joinFirstN : sVector.size();
  if (joinFirstN == 0) {
    return result;
  }
  result = sVector[0];
  for (unsigned int i = 1; i < joinFirstN; ++i) {
    result += "/" + sVector[i];
  }
  return result;
}

/**
 * No longer used in this file and anywhere else in utils. However, it is a public function in a header.
 * Setting it as deprecated instead of removing it.
 */
[[deprecated]] inline std::vector<std::shared_ptr<Scine::Core::Calculator>>
determineUnderlyingCalculators(std::vector<std::string>& listOfMethods, std::vector<std::string>& listOfPrograms) {
  auto& manager = Scine::Core::ModuleManager::getInstance();
  if (listOfMethods.size() != listOfPrograms.size()) {
    throw std::runtime_error("Unequal number of method families and programs given. Please provide a corresponding "
                             "program for each method family in the embedding calculation.");
  }
  if (listOfMethods.size() >= 2 && listOfPrograms.size() >= 2) {
    std::vector<std::shared_ptr<Scine::Core::Calculator>> underlyingCalculators;
    for (long unsigned int i = 0; i < listOfPrograms.size(); i++) {
      underlyingCalculators.emplace_back(manager.get<Scine::Core::Calculator>(
          Scine::Core::Calculator::supports(listOfMethods.at(i)), listOfPrograms.at(i)));
    }
    return underlyingCalculators;
  }
  else {
    throw std::runtime_error(
        "Please provide at least two method families (and the corresponding programs) for an embedding calculation.");
  }
}

inline std::vector<std::shared_ptr<Scine::Core::Calculator>>
getUnderlyingCalculatorsForQMQM(const std::string& multipleMethodsString, const std::string& program) {
  auto& manager = Scine::Core::ModuleManager::getInstance();
  const auto listOfMethods = split(multipleMethodsString, '|');
  std::vector<std::shared_ptr<Scine::Core::Calculator>> underlyingCalculators;
  for (const auto& method : listOfMethods) {
    underlyingCalculators.emplace_back(manager.get<Scine::Core::Calculator>(Scine::Core::Calculator::supports(method), program));
  }
  return underlyingCalculators;
}

inline bool isMMMethod(const std::string& method) {
  std::string small = method;
  boost::algorithm::to_lower(small);
  const std::vector<std::string> mmMethods = {"gaff", "sfam"};
  return std::find(mmMethods.begin(), mmMethods.end(), small) != mmMethods.end();
}

inline bool containsMMMethod(const std::vector<std::string>& methods) {
  for (const auto& method : methods) {
    if (isMMMethod(method)) {
      return true;
    }
  }
  return false;
}

/**
 * @brief Getter for the underlying calculator based on a string representing the model family and a string for the
 * program.
 *
 * Example method families/program combinations:
 * - getCalculator("dft", "orca") - get Orca DFT calculator.
 * - getCalculator("dft", "serenity") - get Serenity DFT calculator.
 * - getCalculator("dft/gaff", "serenity/swoose") - get QM/MM calculator for DFT for the QM region run by Serenity and
 * GAFF/AMBER in the environment, run by Swoose.
 * - getCalculator("cc/dft/gaff", "serenity/serenity/swoose") - get QM/QM/MM embedding.
 * - getCalculator("cc|cc|cc/dft/gaff", "serenity/serenity/swoose") - get QM/QM/MM embedding with CC-in-DFT-in-MM and
 * multiple CC cores.
 *
 * @param method_family The method family string.
 * @param program       The program string.
 * @return The calculator.
 */
inline std::shared_ptr<Scine::Core::Calculator> getCalculator(std::string method_family, std::string program) {
  // Load module manager
  auto& manager = Scine::Core::ModuleManager::getInstance();

  inputPreparation(method_family, program);

  // Generate calculator
  const char sep = '/';
  auto listOfMethods = split(method_family, sep);
  auto listOfPrograms = split(program, sep);

  std::shared_ptr<Scine::Core::Calculator> calc;
  if (listOfMethods.size() < 2 && listOfPrograms.size() < 2) {
    std::shared_ptr<Scine::Core::Calculator> calc;
    try {
      auto f = Scine::Core::Calculator::supports(method_family);
      calc = manager.get<Scine::Core::Calculator>(Scine::Core::Calculator::supports(method_family), program);
    }
    catch (...) {
      if (program.empty() || program == "Any") {
        std::cout << "No SCINE module providing '" << method_family << "' is currently loaded.\n";
        std::cout << "Please add the module to the SCINE_MODULE_PATH\n";
        std::cout << "or load the corresponding Python module in order for it to be accessible.\n";
        throw std::runtime_error("Failed to load method/program.");
      }

      std::cout << "No SCINE module named '" << program << "' providing '" << method_family << "' is currently loaded.\n";
      std::cout << "Please add the module to the SCINE_MODULE_PATH\n";
      std::cout << "or load the corresponding Python module in order for it to be accessible.\n";
      throw std::runtime_error("Failed to load method/program.");
    }
    // Return Calculator
    return calc;
  }
  else {
    if (listOfPrograms.size() != listOfMethods.size()) {
      std::cout << "The number of programs and methods must be identical for embedding calculations" << std::endl;
      throw std::runtime_error("Failed to load method/program.");
    }
    /*
     * Embedding/Multilevel method.
     * Do the following recursion:
     * 1. Determine the type of the embedding calculator from the last method in the method string.
     * 2. Get the embedding calculator QM/MM or QM/QM.
     * 3. Call getCalculator to get the calculator for the environment system QM or MM (these may be multiple in case of
     * QM/QM).
     * 4. Remove the last method from the method string and the last program from the program string.
     * 5. Call getCalculator with the reduced method and program strings to get the underlying calculators (these may be
     * multiple in case of QM/QM).
     */
    const std::string& lastMethod = listOfMethods[listOfMethods.size() - 1];
    const std::string& lastProgram = listOfPrograms[listOfPrograms.size() - 1];
    bool checkForMultipleQMQMRegions = false;
    if (isMMMethod(lastMethod)) {
      calc = manager.get<Scine::Core::Calculator>(Scine::Core::Calculator::supports("QMMM"), "Swoose");
    }
    else if (listOfMethods.size() == 2) {
      // It must be QM/QM or an input error. Try to load Serenity.
      calc = manager.get<Scine::Core::Calculator>(Scine::Core::Calculator::supports("QMQM"), "Serenity");
      checkForMultipleQMQMRegions = true;
    }
    auto castedCalc = std::dynamic_pointer_cast<Scine::Core::EmbeddingCalculator>(calc);
    if (!castedCalc) {
      throw std::runtime_error("Unable to derive the type of embedding calculator. The method family selected is: " + lastMethod);
    }
    const std::string reducedMethodString = join(listOfMethods, listOfMethods.size() - 1);
    const std::string reducedProgramString = join(listOfPrograms, listOfPrograms.size() - 1);
    std::vector<std::shared_ptr<Scine::Core::Calculator>> environmentCalculators;
    std::vector<std::shared_ptr<Scine::Core::Calculator>> activeCalculators;
    if (checkForMultipleQMQMRegions) {
      /*
       * If the QM/QM embedding consists of multiple subsystems, this can be denoted with a '|' to delimit the
       * method families for the individual subsystems. The '/' then separates active and environment subsystems.
       * cc|cc/dft|dft.
       */
      environmentCalculators = getUnderlyingCalculatorsForQMQM(lastMethod, lastProgram);
      activeCalculators = getUnderlyingCalculatorsForQMQM(reducedMethodString, reducedProgramString);
    }
    else {
      environmentCalculators = {getCalculator(lastMethod, lastProgram)};
      activeCalculators = {getCalculator(reducedMethodString, reducedProgramString)};
    }
    std::vector<std::shared_ptr<Scine::Core::Calculator>> listOfCalculators = activeCalculators;
    listOfCalculators.insert(listOfCalculators.end(), environmentCalculators.begin(), environmentCalculators.end());
    castedCalc->setUnderlyingCalculators(listOfCalculators);
    return castedCalc;
  }
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
inline std::shared_ptr<Core::Calculator> spinPropensity(Core::Calculator& calculator, Core::Log& log, int maxSpinDifference) {
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
    clone->settings().modifyInt(Utils::SettingsNames::spinMultiplicity, i);
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
                << std::to_string(bestClone->settings().getInt(Utils::SettingsNames::spinMultiplicity)) << Core::Log::nl;
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
  // clang-format off
  std::vector<std::string> exceptions = {
      "PNO-CC",
      "HF-3C",
      "PBEH-3C",
      "B97-3C",
      "-F12",
  };
  // if input contains one of these as a substring
  // we ignore the number of '-' in the exception and still check for dispersion
  std::vector<std::string> exceptionsWithPotentialDisp = {
      "CAM-B3LYP",
      "M05-2X",
      "M06-L",
      "M06-2X",
      "M06-HF",
      "M08-HX",
      "M08-SO",
      "M11-L",
      "MN12-L",
      "MN12-SX",
      "MN15-L",
      "LC-PBE",
      "LC-WPBE",
  };
  // clang-format on
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
  for (const auto& exceptionMethod : exceptionsWithPotentialDisp) {
    if (inputCopy.find(exceptionMethod) != std::string::npos) {
      auto allowed_hyphens = std::count(exceptionMethod.begin(), exceptionMethod.end(), '-');
      std::vector<std::string> newSegments;
      newSegments.push_back("");
      for (int i = 0; i <= allowed_hyphens; ++i) {
        if (i == 0)
          newSegments[0] = segments[i];
        else
          newSegments[0] += "-" + segments[i];
      }
      for (unsigned long i = allowed_hyphens + 1; i < segments.size(); ++i) {
        newSegments.push_back(segments[i]);
      }
      segments = newSegments;
      break;
    }
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

/*
 * @brief Returns results for zero electron systems if energy should be zero (no empirics)
 * @param atoms: the atoms (determines the sizes of results and charges)
 * @param requiredProperties: the properties that should be filled
 */
inline Results calculateZeroElectrons(const AtomCollection& atoms, const PropertyList& requiredProperties) {
  int nAtoms = atoms.size();
  if (nAtoms > 1) {
    // we deem multi-atom-no-electron systems not important enough
    // to calculate nuclear repulsion energy + gradients ourselves
    throw std::runtime_error("Calculation of multinuclear system with zero electrons is not implemented.");
  }
  Results results;
  if (requiredProperties.containsSubSet(Property::Energy)) {
    results.set<Property::Energy>(0.0);
  }
  if (requiredProperties.containsSubSet(Property::Gradients)) {
    results.set<Property::Gradients>(GradientCollection::Zero(nAtoms, 3));
  }
  if (requiredProperties.containsSubSet(Property::Hessian)) {
    results.set<Property::Hessian>(HessianMatrix(nAtoms, nAtoms));
  }
  if (requiredProperties.containsSubSet(Property::BondOrderMatrix)) {
    results.set<Property::BondOrderMatrix>(BondOrderCollection(nAtoms));
  }
  if (requiredProperties.containsSubSet(Property::AtomicCharges)) {
    std::vector<double> charges;
    for (const auto& e : atoms.getElements()) {
      charges.push_back(ElementInfo::Z(e));
    }
    results.set<Property::AtomicCharges>(charges);
  }
  results.set<Property::SuccessfulCalculation>(true);
  return results;
}

} // namespace CalculationRoutines
} // namespace Utils
} // namespace Scine

#endif // UTILS_CALCULATION_ROUTINES_H_
