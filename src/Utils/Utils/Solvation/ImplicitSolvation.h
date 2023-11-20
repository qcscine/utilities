/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_IMPLICIT_SOLVATION_H_
#define UTILS_IMPLICIT_SOLVATION_H_

#include "Utils/Settings.h"
#include "Utils/UniversalSettings/SettingsNames.h"
#include <iostream>
#include <string>

namespace Scine {
namespace Utils {
namespace Solvation {

namespace ImplicitSolvation {

/**
 * @brief Checks solvation and solvent in settings to be a valid input with the given available models
 *
 * @note if solvation == 'any': the first available solvation model is applied
 * @note if solvent == 'any': 'water' is applied and a warning is given in std::out
 *
 * @param availableSolvationModels List of available solvation models of the calculator
 * @param settings The calculator settings to get the solvation and solvent input, 'any' entries are updated
 *
 * @throws std::logic_error for wrong input
 * @return bool whether implicit solvation has to be applied or not
 */
inline bool solvationNeededAndPossible(std::vector<std::string> availableSolvationModels, Settings& settings) {
  std::string solvent = settings.getString(Utils::SettingsNames::solvent);
  std::string solvation = settings.getString(Utils::SettingsNames::solvation);
  /* make all strings lowercase */
  std::for_each(solvent.begin(), solvent.end(), [](char& c) { c = ::tolower(c); });
  std::for_each(solvation.begin(), solvation.end(), [](char& c) { c = ::tolower(c); });
  for (std::string& model : availableSolvationModels) {
    std::for_each(model.begin(), model.end(), [](char& c) { c = ::tolower(c); });
  }

  // no solvation
  if (solvation == "none" || solvation.empty()) {
    // but solvent has been set -> error
    if (solvent != "none" && !solvent.empty()) {
      throw std::logic_error("Error: Specified no solvation, but specified solvent " + solvent);
    }
    // no solvation, no solvent -> valid input and no solvation needed
    return false;
  }
  // solvation is not None, and calculator does not have any available solvation models -> error
  if (availableSolvationModels.empty()) {
    throw std::logic_error("Error: You specified the solvation model '" + solvation +
                           "'. However, this calculator does not support solvation.");
  }
  // calculator supports solvation, but solvation is something that is not supported -> error
  if (solvation != "any" && std::find(availableSolvationModels.begin(), availableSolvationModels.end(), solvation) ==
                                availableSolvationModels.end()) {
    /* prepare output */
    std::string modelString;
    for (const std::string& model : availableSolvationModels) {
      modelString += model + "\n";
    }
    throw std::logic_error("Error: You specified the solvation model '" + solvation +
                           "'. This is not possible with this calculator.\n"
                           "The available solvation models are:\n" +
                           modelString);
  }
  // solvation is valid
  // but solvent is specified to be none -> error
  if (solvent.empty() || solvent == "none") {
    throw std::logic_error("Error: Specified solvation model " + solvation + " but specified solvent as none");
  }
  // solvent is specified as 'any' -> warning and select water
  if (solvent == "any") {
    // TODO proper warning
    std::cerr << "Warning, specified implicit solvation with '" + solvation + "', but solvent was set to 'any'. Using water as default."
              << std::endl;
    solvent = "water";
  }
  settings.modifyString(Utils::SettingsNames::solvent, solvent);
  // updating 'any' entry with actual model
  if (solvation == "any") {
    solvation = availableSolvationModels.front();
  }
  settings.modifyString(Utils::SettingsNames::solvation, solvation);
  // all harmful things should have been caught -> enabling solvation
  return true;
}

} /* namespace ImplicitSolvation */
} /* namespace Solvation */
} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_IMPLICIT_SOLVATION_H_
