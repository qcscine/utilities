/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Geometry/FormulaGenerator.h"
#include "Utils/Geometry/ElementInfo.h"
#include <map>

namespace Scine {
namespace Utils {

std::string generateChemicalFormula(const ElementTypeCollection& elements, const std::string& numberPrefix,
                                    const std::string& numberPostfix) {
  if (elements.empty())
    return "(empty)";

  // Treat C and H separately
  unsigned nC = 0;
  unsigned nH = 0;

  std::map<std::string, unsigned> elementsMap;
  for (auto e : elements) {
    if (e == ElementType::C)
      ++nC;
    else if (e == ElementType::H)
      ++nH;
    else {
      std::string elementSymbol = ElementInfo::symbol(e);
      elementsMap[elementSymbol]++;
    }
  }

  std::string string;
  if (nC != 0)
    string += singleElementPartOfFormula("C", nC, numberPrefix, numberPostfix);
  if (nH != 0)
    string += singleElementPartOfFormula("H", nH, numberPrefix, numberPostfix);
  for (const auto& s : elementsMap) {
    string += singleElementPartOfFormula(s.first, s.second, numberPrefix, numberPostfix);
  }

  return string;
}

std::string singleElementPartOfFormula(std::string symbol, int numberOccurrences, const std::string& numberPrefix,
                                       const std::string& numberPostfix) {
  if (numberOccurrences != 1) {
    symbol += numberPrefix;
    symbol += std::to_string(numberOccurrences);
    symbol += numberPostfix;
  }
  return symbol;
}

} /* namespace Utils */
} /* namespace Scine */
