/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 * @brief This header file contains functions that convert structure data into chemical formulas.
 */
#ifndef UTILS_FORMULAGENERATOR_H_
#define UTILS_FORMULAGENERATOR_H_

#include "Utils/Typenames.h"
#include <string>

namespace Scine {
namespace Utils {

/**
 * @brief This function returns the elemental composition of a compound based on its ElementTypeCollection.
 *
 * The convention is the Hill Order System.
 * It starts with C, then H, then all the other elements in alphabetical order.
 * Returns "(empty)" if elements is empty.
 *
 * @param elements The collection of elements to be transcribed.
 * @param numberPrefix Characters to put before each number for rich text formats (Latex etc).
 * @param numberPostfix Characters to put after each number for rich text formats (Latex etc).
 * @return std::string The final string.
 */
std::string generateChemicalFormula(const ElementTypeCollection& elements, const std::string& numberPrefix = "",
                                    const std::string& numberPostfix = "");
/**
 * @brief This function transcribes a single element into a string.
 *
 * @param symbol The element symbol
 * @param numberOccurrences The number of times the element is present in a given structure.
 * @param numberPrefix Characters to put before each number for rich text formats (Latex etc).
 * @param numberPostfix Characters to put after each number for rich text formats (Latex etc).
 * @return std::string The final string.
 */
std::string singleElementPartOfFormula(std::string symbol, int numberOccurrences, const std::string& numberPrefix,
                                       const std::string& numberPostfix);

} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_FORMULAGENERATOR_H_
