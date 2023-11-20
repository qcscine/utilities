/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UNIVERSALSETTINGS_DESCRIPTORCOLLECTION_H
#define UNIVERSALSETTINGS_DESCRIPTORCOLLECTION_H
/* Internal Headers */
#include "Utils/UniversalSettings/GenericDescriptor.h"
#include "Utils/UniversalSettings/SettingDescriptor.h"
/* External Headers */
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace Scine {
namespace Utils {
namespace UniversalSettings {
/* Forward Declarations */
class ValueCollection;

/*!
 * @brief Setting that wraps a vector<pair<str, GenericDescriptor>>.
 *
 * A setting (string descriptor, generic value pair) whose value is a list of
 * string, GenericDescriptor pairs, which altogether define a nested list of
 * descriptors.
 */
class DescriptorCollection : public SettingDescriptor {
 public:
  DescriptorCollection() = default;
  explicit DescriptorCollection(std::string description);

  std::unique_ptr<SettingDescriptor> clone() const override;
  GenericValue getDefaultGenericValue() const override;
  bool validValue(const GenericValue& v) const override;
  std::string explainInvalidValue(const GenericValue& v) const override;
  std::string explainInvalidValue(const ValueCollection& v) const;
  std::map<std::string, std::string> gatherInvalidExplanations(const ValueCollection& v) const;
  std::map<std::string, std::string> gatherInvalidExplanations(const GenericValue& v) const;

  using KeyValuePair = std::pair<std::string, GenericDescriptor>;
  using Container = std::vector<KeyValuePair>;
  using value_type = Container::value_type;
  using iterator = Container::iterator;
  using const_iterator = Container::const_iterator;

  //!@name Iterators
  //!@{
  iterator begin();
  const_iterator begin() const;
  iterator end();
  const_iterator end() const;
  //!@}

  //!@name Accessors and modifiers
  //!@{
  void push_back(std::string key, GenericDescriptor e);
  GenericDescriptor& operator[](unsigned i);
  const GenericDescriptor& operator[](unsigned i) const;
  GenericDescriptor& at(unsigned i);
  const GenericDescriptor& at(unsigned i) const;
  GenericDescriptor& get(const std::string& key);
  const GenericDescriptor& get(const std::string& key) const;
  //!@}

  //!@name Information
  //!@{
  //! Checks whether a particular descriptor field exists in the configuration
  bool exists(const std::string& key) const;
  //! Checks whether the list of string, descriptor pairs is empty
  bool empty() const;
  //! Returns the number of string, descriptor pairs
  int size() const;

  /**
   * @brief Determines whether a corresponding list of string, GenericValue
   *   pairs matches the configuration of this setting.
   *
   * @param v A list of string, GenericValue pairs
   *
   * @return Whether the supplied values match the setting configuration
   */
  bool validValue(const ValueCollection& v) const;
  //!@}

 private:
  Container descriptors_;
};

std::string invalidSettingsMapToString(const std::map<std::string, std::string>& invalidSettings);
ValueCollection createDefaultValueCollection(const DescriptorCollection& descriptors);

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */

#endif // UNIVERSALSETTINGS_DESCRIPTORCOLLECTION_H
