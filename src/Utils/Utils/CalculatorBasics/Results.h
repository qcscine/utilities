/**
 * @file
 * @brief A file containing definitions of the results
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_RESULTS_H
#define UTILS_RESULTS_H

#include "PropertyList.h"
#include <boost/any.hpp>
#include <exception>
#include <memory>
#include <string>

namespace Scine {
namespace Utils {

/**
 * Exception thrown if a property is requested and it is not stored.
 */
class PropertyNotPresentException : public std::exception {
 public:
  const char* what() const noexcept final {
    return "Property desired not present in results.";
  }
};

/**
 * Class for the properties obtained in a single-point calculation.
 * To obtain the properties:
 * - "get" methods return references to the results (also allows simple copy)
 * - "take" methods move the results (which will not be present in Results afterwards).
 *
 * Template pattern does not allow for pImpl idiom without defining a specialized template for each Property.
 */
class Results {
 public:
  Results();
  ~Results();
  Results(Results&& rhs) noexcept;
  Results& operator=(Results&& rhs) noexcept;
  Results(const Results& rhs);
  Results& operator=(const Results& rhs);

  template<Property property>
  bool has() const {
    return resultsMap_.find(property) != resultsMap_.end();
  };

  /** \brief Function that returns PropertyList of all Properties contained in the results */
  PropertyList allContainedProperties() const;

  template<Property property>
  const typename PropertyType<property>::Type& get() const {
    if (!has<property>()) {
      throw PropertyNotPresentException();
    }
    // need to cast boost any to right type
    return boost::any_cast<const typename PropertyType<property>::Type&>(resultsMap_.at(property));
  }
  template<Property property>
  void set(typename PropertyType<property>::Type dataToSet) {
    resultsMap_[property] = dataToSet;
  }
  template<Property property>
  typename PropertyType<property>::Type take() {
    if (!has<property>()) {
      throw PropertyNotPresentException();
    }
    auto result = boost::any_cast<typename PropertyType<property>::Type>(resultsMap_.at(property));
    resultsMap_.erase(property);
    return result;
  }

 private:
  std::map<Property, boost::any> resultsMap_;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_RESULTS_H
