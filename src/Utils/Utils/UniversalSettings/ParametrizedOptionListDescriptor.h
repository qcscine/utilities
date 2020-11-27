/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UNIVERSALSETTINGS_PARAMETRIZEDOPTIONLISTDESCRIPTOR_H
#define UNIVERSALSETTINGS_PARAMETRIZEDOPTIONLISTDESCRIPTOR_H
/* Internal Headers */
#include "Utils/UniversalSettings/DescriptorCollection.h"
#include "Utils/UniversalSettings/SettingDescriptor.h"
/* External Headers */
#include <string>
#include <utility>
#include <vector>

namespace Scine {
namespace Utils {
namespace UniversalSettings {

/*!
 * @brief SettingDescriptor for a list of options, of which one must be chosen,
 *   with corresponding settings (that depend on the exact option).
 *
 * This is for example useful if one of several algorithms must be chosen for
 * some task and each algorithm has specific settings.
 */
class ParametrizedOptionListDescriptor : public SettingDescriptor {
 public:
  using OptionAndSettings = std::pair<std::string, DescriptorCollection>;

  /**
   * @brief Constructor
   *
   * @param propertyDescription A string describing what the options are for
   */
  explicit ParametrizedOptionListDescriptor(std::string propertyDescription);

  std::unique_ptr<SettingDescriptor> clone() const override;
  GenericValue getDefaultGenericValue() const override;
  bool validValue(const GenericValue& v) const override;

  //! Returns the number of available options
  int optionCount() const;
  //! Returns whether a string-identified option exists
  bool optionExists(const std::string& name) const;
  /*!
   * @brief Returns the index of the default option within the list of options
   *
   * @throws EmptyOptionListException if the list of options is empty
   */
  int getDefaultIndex() const;
  //! Returns the default chosen value
  const std::string& getDefaultValue() const;
  //! Returns the default chosen value
  const std::string& getDefaultOption() const;
  //! Read-only accessor to list of option strings and their matching settings
  const std::vector<OptionAndSettings>& getAllOptions() const;

  /*! Add an option for which there are no specific settings. */
  void addOption(std::string option);
  /*! Add an option with the corresponding specific settings. */
  void addOption(std::string option, DescriptorCollection optionSettings);
  /**
   * @brief Sets the default option to a particular string
   *
   * @param def The option to select as default
   *
   * @throws OptionDoesNotExistException if the supplied option does not match
   *   an existing option
   */
  void setDefaultOption(const std::string& def);

  //! Fetch the settings corresponding to a particular option
  const DescriptorCollection& getSettings(const std::string& option) const;
  //! Fetch the settings corresponding to the default option
  const DescriptorCollection& getDefaultSettings() const;

 private:
  int getIndex(const std::string& option) const;
  std::vector<OptionAndSettings> options_;
  int defaultIndex_ = 0;
};

inline const std::string& ParametrizedOptionListDescriptor::getDefaultValue() const {
  return getDefaultOption();
}

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */
#endif // UNIVERSALSETTINGS_PARAMETRIZEDOPTIONLISTDESCRIPTOR_H
