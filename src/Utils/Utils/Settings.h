/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
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
 *  - A derived class of the Settings class has to define a construcor
 *    that populates the protected _fields member with a set of possible
 *    options and their default/descriptions.
 */
class Settings : public Scine::Utils::UniversalSettings::ValueCollection {
 public:
  /**
   * @brief Construct a new Settings object
   *
   * @param name
   */
  explicit Settings(std::string name) : _name(name), _fields(Utils::UniversalSettings::DescriptorCollection(name)) {
  }
  /**
   * @brief Constructor from ValueCollection and DescriptorCollection.
   */
  Settings(Utils::UniversalSettings::ValueCollection settings, Utils::UniversalSettings::DescriptorCollection fields)
    : Utils::UniversalSettings::ValueCollection(std::move(settings)), _fields(std::move(fields)) {
  }
  /**
   * @brief Getter for the name set in the constructor of the settings object.
   * @return std::string Returns the name.
   */
  const std::string& name() const {
    return _name;
  }
  /**
   * @brief Checks if the current settings are acceptable w.r.t. the
   *        defined boundaries.
   *
   * @return true  If the settings are valid.
   * @return false If the settings are not valid.
   */
  bool check() const {
    const Scine::Utils::UniversalSettings::ValueCollection& self = *this;
    return _fields.validValue(self);
  }

  /**
   * @brief Resets the current settings to the default ones.
   */
  void resetToDefaults() {
    Scine::Utils::UniversalSettings::ValueCollection& self = *this;
    self = Utils::UniversalSettings::createDefaultValueCollection(this->_fields);
  }

  /**
   * @brief Returns a const reference of the protected member _fields.
   * @return const Utils::UniversalSettings::DescriptorCollection& Underlying descriptor collection.
   */
  const Utils::UniversalSettings::DescriptorCollection& getDescriptorCollection() const {
    return _fields;
  }

  /**
   * @brief Returns the default ValueCollection as a settings object for a given key and option.
   * @param key The key describing the ParametrizedOptionListDescriptor.
   * @param option The option of which the default ValueCollection should be returned.
   * @return Settings The default settings corresponding to the option.
   */
  Settings getDefaultSettingsForOptionListWithSettings(const std::string& key, const std::string& option) {
    if (_fields.exists(key)) {
      if (_fields.get(key).relatesToParametrizedOptionList()) {
        auto optionList = _fields.get(key).getParametrizedOptionListDescriptor();
        for (const auto& o : optionList.getAllOptions()) {
          if (o.first == option)
            return Settings(Utils::UniversalSettings::createDefaultValueCollection(o.second), o.second);
        }
        throw std::runtime_error("The given option is not a valid option.");
      }
      else {
        throw std::runtime_error("The given setting is not an option list with settings.");
      }
    }
    else {
      throw Core::SettingsKeyError();
    }
  }

 protected:
  /**
   * @brief The blueprint for the allowed options.
   *
   * All derived classes need to populate this collection
   * in their constructor.
   */
  std::string _name;
  Utils::UniversalSettings::DescriptorCollection _fields;
};

} // namespace Utils
} // namespace Scine

#endif
