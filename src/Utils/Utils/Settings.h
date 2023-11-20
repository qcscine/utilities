/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_SETTINGS_H_
#define UTILS_SETTINGS_H_

#include <Core/Exceptions.h>
#include <Utils/UniversalSettings/DescriptorCollection.h>
#include <Utils/UniversalSettings/ParametrizedOptionListDescriptor.h>
#include <Utils/UniversalSettings/ValueCollection.h>

namespace Scine {
namespace Utils {

/**
 * @class Settings Settings.h
 * @brief An interface for settings of all sorts.
 *
 * This interface allows all other interfaces in Scine::Core to
 * have a central place for Settings, while their implementations
 * maintain their specific settings independently.
 *
 * Implementation Notes:
 *  - A derived class of the Settings class has to define a constructor
 *    that populates the protected _fields member with a set of possible
 *    options and their default/descriptions.
 */
class Settings : public UniversalSettings::ValueCollection {
 public:
  /*! Normalize a string values' case if a case-insensitive option can be
   *   matched in the corresponding descriptor
   *
   * For OptionList- and ParametrizedOptionListDescriptors, if the current
   * option value identifier does not exactly match any option within the
   * descriptor, a case-insensitive match is sought. If one is found, the
   * value identifier is overwritten with the matched descriptor identifier.
   *
   * @note Recursively descends into all nested
   *   DescriptorCollection/ValueCollection pairs.
   */
  static void normalizeStringCases(const UniversalSettings::DescriptorCollection& descriptors, ValueCollection& values);

  //!@name Constructors
  //!@{
  Settings() = default;
  /**
   * @brief Construct a new Settings object
   *
   * @param name
   */
  explicit Settings(std::string name) : _name(std::move(name)), _fields(UniversalSettings::DescriptorCollection(name)) {
  }
  /**
   * @brief Constructor from ValueCollection and DescriptorCollection.
   */
  Settings(UniversalSettings::ValueCollection settings, UniversalSettings::DescriptorCollection fields)
    : UniversalSettings::ValueCollection(std::move(settings)),
      _name(fields.getPropertyDescription()),
      _fields(std::move(fields)) {
  }
  //!@}

  /**
   * @brief Virtual destructor.
   */
  virtual ~Settings() = default;

  //!@name Modification
  //!@{
  /**
   * @brief Merges values of matching keys from a collection into settings
   *
   * @param collection Values to merge if keys match
   * @param allowSuperfluous If true, the function will ignore value keys that are not present in the descriptors.
   *                         If false, the function will throw an error if such a case is encountered.
   */
  void merge(const Utils::UniversalSettings::ValueCollection& collection, bool allowSuperfluous = false) {
    for (const auto& [key, value] : collection) {
      if (valueExists(key)) {
        modifyValue(key, value);
      }
      else if (!allowSuperfluous) {
        throw std::logic_error("Error: the key: '" + key + "' was not recognized.");
      }
    }
  }
  //!@name Modification
  //!@{
  /**
   * @brief Merges two collections, if key of given collection does not exist, it is added.
   *
   * @param collection Values to merge
   */
  void mergeAll(const Utils::UniversalSettings::ValueCollection& collection) {
    for (const auto& [key, value] : collection) {
      if (valueExists(key)) {
        modifyValue(key, value);
      }
      else {
        addGenericValue(key, value);
      }
    }
  }

  //!@overload
  void normalizeStringCases();

  /**
   * @brief Resets the current settings to the default ones.
   */
  void resetToDefaults() {
    UniversalSettings::ValueCollection& self = *this;
    self = UniversalSettings::createDefaultValueCollection(this->_fields);
  }
  //!@}

  //!@name Information
  //!@{
  /**
   * @brief Getter for the name set in the constructor of the settings object.
   * @return std::string Returns the name.
   */
  const std::string& name() const {
    return _name;
  }

  [[deprecated("Prefer API-equivalent 'valid' function")]] bool check() const {
    return valid();
  }

  /**
   * @brief Checks if the current settings are acceptable w.r.t. the
   *        defined boundaries.
   *
   * @returns Whether the settings are valid w.r.t their descriptor bounds
   */
  bool valid() const {
    return _fields.validValue(*this);
  }

  /**
   * @brief Throw an exception containing explanations for all setting values
   *   that violate their descriptor bounds.
   *
   * @pre A setting value is invalid w.r.t. descriptor bounds, i.e.
   *   `this->valid() == false`
   */
  void throwIncorrectSettings() const {
    assert(!_fields.validValue(*this));
    throw UniversalSettings::InvalidSettingsException(_fields.explainInvalidValue(*this));
  }

  /**
   * @brief Returns a const reference of the protected member _fields.
   * @return const UniversalSettings::DescriptorCollection& Underlying descriptor collection.
   */
  const UniversalSettings::DescriptorCollection& getDescriptorCollection() const {
    return _fields;
  }

  /**
   * @brief Returns the default ValueCollection as a settings object for a given key and option.
   * @param key The key describing the ParametrizedOptionListDescriptor.
   * @param option The option of which the default ValueCollection should be returned.
   * @return Settings The default settings corresponding to the option.
   */
  Settings getDefaultSettingsForOptionListWithSettings(const std::string& key, const std::string& option) const {
    if (_fields.exists(key)) {
      if (_fields.get(key).relatesToParametrizedOptionList()) {
        auto optionList = _fields.get(key).getParametrizedOptionListDescriptor();
        for (const auto& o : optionList.getAllOptions()) {
          if (o.first == option) {
            return Settings(UniversalSettings::createDefaultValueCollection(o.second), o.second);
          }
        }
        throw std::runtime_error("The given option is not a valid option.");
      }

      throw std::runtime_error("The given setting is not an option list with settings.");
    }

    throw Core::SettingsKeyError();
  }
  //!@}

 protected:
  //! Name of the settings object
  std::string _name;
  /**
   * @brief The blueprint for the allowed options.
   *
   * All derived classes need to populate this collection
   * in their constructor.
   */
  UniversalSettings::DescriptorCollection _fields;
};

} // namespace Utils
} // namespace Scine

#endif
