/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/CalculatorBasics/Results.h"
#include "PropertyList.h"
#include "Utils/Bonds/BondOrderCollection.h"
#include <boost/optional.hpp>

namespace Scine {
namespace Utils {

Results::Results() = default;
Results::Results(Results&& /*rhs*/) noexcept = default;
Results::Results(const Results& rhs) = default;

Results& Results::operator=(const Results& rhs) {
  static_assert(std::is_move_assignable<Results>::value, "Must be move assignable!");
  *this = Results(rhs);
  return *this;
}

Results Results::operator+(const Results& rhs) const {
  auto combined = Results(*this);
  combined += rhs;
  return combined;
}

Results& Results::operator+=(const Results& rhs) {
  // Loop manually, because we are lacking definitions for merge
  for (auto property : allProperties) {
    const auto rhsIt = rhs.resultsMap_.find(property);
    if (rhsIt != rhs.resultsMap_.end()) {
      this->resultsMap_[rhsIt->first] = rhsIt->second;
    }
  }
  return *this;
}

PropertyList Results::allContainedProperties() const {
  PropertyList propertyList;
  // Loop over all possible properties and check whether they are available
  for (auto property : allProperties) {
    if (resultsMap_.find(property) != resultsMap_.end()) {
      propertyList.addProperty(property);
    }
  }
  return propertyList;
}

Results::~Results() = default;

Results& Results::operator=(Results&& /*rhs*/) noexcept = default;

} // namespace Utils
} // namespace Scine
