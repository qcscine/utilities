/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Technical/ScopedLocale.h"
#include <locale>

namespace Scine {
namespace Utils {

ScopedLocale::ScopedLocale(const std::string& targetLocale) {
  originalLocale_ = std::locale("").name();
  std::locale::global(std::locale(targetLocale.c_str()));
}

ScopedLocale::~ScopedLocale() {
  if (mustBeSetBackInDestructor_) {
    std::locale::global(std::locale(originalLocale_.c_str()));
  }
}

ScopedLocale::ScopedLocale(ScopedLocale&& other) noexcept {
  originalLocale_ = std::move(other.originalLocale_);
  other.mustBeSetBackInDestructor_ = false;
}

ScopedLocale& ScopedLocale::operator=(ScopedLocale&& other) noexcept {
  originalLocale_ = std::move(other.originalLocale_);
  other.mustBeSetBackInDestructor_ = false;
  mustBeSetBackInDestructor_ = true;
  return *this;
}

ScopedLocale ScopedLocale::cLocale() {
  return ScopedLocale{"C"};
}

} /* namespace Utils */
} /* namespace Scine */
