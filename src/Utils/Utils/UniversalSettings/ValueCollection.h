/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UNIVERSALSETTINGS_VALUECOLLECTION_H
#define UNIVERSALSETTINGS_VALUECOLLECTION_H

/* Internal Headers */
#include "Utils/UniversalSettings/Exceptions.h"
#include "Utils/UniversalSettings/GenericValue.h"
#include "Utils/UniversalSettings/OptionListDescriptor.h"
/* External Headers */
#include <string>
#include <utility>
#include <vector>

namespace Scine {
namespace Utils {
namespace UniversalSettings {

/*!
 * @brief Wrapper around vector<pair<string, GenericValue>>
 *
 * Class holding values corresponding to some SettingDescriptorCollection.
 * The functions addXXX() throw a AlreadyExistingValueException if the given name already exists.
 * The functions getXXX() throw a InexistingValueException if the given value does not exist, and a
 * InvalidValueConversionException if the incorrect type is requested.
 * The functions modifyXXX() assume the key is already existing and throw a InexistingValueException if not,
 * and a InvalidValueConversionException if the type is not the same.
 */
class ValueCollection {
 public:
  using KeyValuePair = std::pair<std::string, GenericValue>;
  using Container = std::vector<KeyValuePair>;

  virtual ~ValueCollection() = default;

  [[deprecated("Prefer equivalent function size")]] int count() const;
  int size() const;
  bool empty() const;

  std::vector<std::string> getKeys() const;
  Container items() const;
  bool valueExists(const std::string& name) const;

  void addGenericValue(std::string name, GenericValue value);
  GenericValue getValue(const std::string& name) const;
  void modifyValue(const std::string& name, GenericValue value);
  void dropValue(const std::string& name);

  /**
   * @brief Extracts a value from the ValueCollection, substitutes a default if not present
   *
   * @note Removes the key-value pair from the ValueCollection if present
   *
   * @tparam T Default value type if the key is not present in the value collection
   * @param name The key to the value to be extracted
   * @param defaultValue Default value to substitute if no key is present
   *
   * @return The value for the key in the collection if present, the default value otherwise
   */
  template<typename T>
  T extract(const std::string& name, T defaultValue) {
    if (valueExists(name)) {
      try {
        defaultValue = getValue(name);
      }
      catch (const std::runtime_error&) {
        throw ValueHasDifferentTypeException(name);
      }
      dropValue(name);
    }
    return defaultValue;
  }

  void addBool(std::string name, bool value);
  bool getBool(const std::string& name) const;
  void modifyBool(const std::string& name, bool value);

  void addInt(std::string name, int value);
  int getInt(const std::string& name) const;
  void modifyInt(const std::string& name, int value);

  void addDouble(std::string name, double value);
  double getDouble(const std::string& name) const;
  void modifyDouble(const std::string& name, double value);

  void addString(std::string name, std::string value);
  std::string getString(const std::string& name) const;
  void modifyString(const std::string& name, std::string value);

  void addCollection(std::string name, ValueCollection value);
  ValueCollection getCollection(const std::string& name) const;
  void modifyCollection(const std::string& name, ValueCollection value);

  void addOptionWithSettings(std::string name, ParametrizedOptionValue value);
  ParametrizedOptionValue getOptionWithSettings(const std::string& name) const;
  void modifyOptionsWithSettings(const std::string& name, ParametrizedOptionValue value);

  void addIntList(std::string name, GenericValue::IntList value);
  GenericValue::IntList getIntList(const std::string& name) const;
  void modifyIntList(const std::string& name, GenericValue::IntList value);

  void addIntListList(std::string name, GenericValue::IntListList value);
  GenericValue::IntListList getIntListList(const std::string& name) const;
  void modifyIntListList(const std::string& name, GenericValue::IntListList value);

  void addDoubleList(std::string name, GenericValue::DoubleList value);
  GenericValue::DoubleList getDoubleList(const std::string& name) const;
  void modifyDoubleList(const std::string& name, GenericValue::DoubleList value);

  void addStringList(std::string name, GenericValue::StringList value);
  GenericValue::StringList getStringList(const std::string& name) const;
  void modifyStringList(const std::string& name, GenericValue::StringList value);

  void addCollectionList(std::string name, GenericValue::CollectionList value);
  GenericValue::CollectionList getCollectionList(const std::string& name) const;
  void modifyCollectionList(const std::string& name, GenericValue::CollectionList value);

  Container::const_iterator begin() const {
    return std::begin(values_);
  }

  Container::const_iterator end() const {
    return std::end(values_);
  }

  std::vector<std::string> getDivergingKeys(const ValueCollection& otherSettings, bool onlyFirstDivergingKey = false) const;

 private:
  const GenericValue& getGenericValue(const std::string& name) const;
  Container::iterator findName(const std::string& name);
  Container::const_iterator findName(const std::string& name) const;

  Container values_;
};

bool operator==(const ValueCollection& lhs, const ValueCollection& rhs);
bool operator!=(const ValueCollection& lhs, const ValueCollection& rhs);

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */

#endif // UNIVERSALSETTINGS_VALUECOLLECTION_H
