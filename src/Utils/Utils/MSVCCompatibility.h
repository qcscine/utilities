/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_MSVCCOMPATIBILITY_H
#define UTILS_MSVCCOMPATIBILITY_H

#include <iso646.h>
#include <array>
#include <fstream>
#include <numeric>

#define _USE_MATH_DEFINES
#include <cmath>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

inline void SetEnv_(std::string p_Tag, std::string p_Value) {
#if defined(_WIN32) || defined(_WIN64)
  int ret = _putenv((p_Tag + p_Value).c_str());
#elif defined(__linux__)
  setenv(p_Tag.c_str(), p_Value.c_str(), 1);
#endif
}

inline void UnsetEnv_(std::string p_Tag) {
#if defined(_WIN32) || defined(_WIN64)
  int ret = _putenv((p_Tag + "=").c_str());
#elif defined(__linux__)
  unsetenv(p_Tag.c_str());
#endif
}

#endif // UTILS_MSVCCOMPATIBILITY_H
