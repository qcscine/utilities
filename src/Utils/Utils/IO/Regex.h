/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_REGEX_H_
#define UTILS_REGEX_H_

#include <string>

namespace Scine {
namespace Utils {

/**
 * @brief Helpers for std::regex.
 *
 * The functions deliver regular expression building blocks to be used when composing regular expressions.
 * Functions starting with "capturing" add a parenthesis to capture the value.
 */
namespace Regex {
/**
 * @brief Regex code for the beginning of a line.
 * @return std::string Returns Regex code for the beginning of a line.
 */
std::string lineBegin();
/**
 * @brief Regex code for the end of a line.
 * @return std::string Returns Regex code for the end of a line.
 */
std::string lineEnd();
/**
 * @brief Regex code for a floating point number.
 * @return std::string Returns Regex code for a floating point number.
 */
std::string floatingPointNumber();
/**
 * @brief Regex code for a floating point number including capture.
 * @return std::string Returns Regex code for a floating point number including capture.
 */
std::string capturingFloatingPointNumber();
/**
 * @brief Regex code for an integer number.
 * @return std::string Returns Regex code for an integer number.
 */
std::string integerNumber();
/**
 * @brief Regex code for an integer number including capture.
 * @return std::string Returns Regex code for an integer number including capture.
 */
std::string capturingIntegerNumber();
/**
 * @brief Regex code for an element symbol.
 * @return std::string Returns Regex code for an element symbol.
 */
std::string elementSymbol();
/**
 * @brief Regex code for an element symbol including capture.
 * @return std::string Returns Regex code for an element symbol including capture.
 */
std::string capturingElementSymbol();
/**
 * @brief Adds capture parentheses to any string.
 * @param s The Regex string describing the type to be captured.
 * @return std::string Returns the string s with capture parenthesis added.
 */
std::string addCaptureParenthesis(const std::string& s);

} /* namespace Regex */
} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_REGEX_H_
