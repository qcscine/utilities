/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/IO/Regex.h"

namespace Scine {
namespace Utils {
namespace Regex {

std::string lineBegin() {
  return R"(^)";
}

std::string lineEnd() {
  return R"($)";
}

std::string floatingPointNumber() {
  // NB: syntax (?:.....), instead of (.....), makes sure that the parenthesis is not captured
  return R"([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)";
}

std::string capturingFloatingPointNumber() {
  return addCaptureParenthesis(floatingPointNumber());
}

std::string integerNumber() {
  return R"([+-]?\d+)";
}

std::string capturingIntegerNumber() {
  return addCaptureParenthesis(integerNumber());
}

std::string addCaptureParenthesis(const std::string& s) {
  return "(" + s + ")";
}

} /* namespace Regex */
} /* namespace Utils */
} /* namespace Scine */
