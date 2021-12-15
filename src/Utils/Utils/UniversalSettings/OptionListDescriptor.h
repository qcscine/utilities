/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UNIVERSALSETTINGS_OPTIONLISTDESCRIPTOR_H
#define UNIVERSALSETTINGS_OPTIONLISTDESCRIPTOR_H
/* Internal Headers */
#include "Utils/UniversalSettings/SettingDescriptor.h"
/* External Headers */
#include <string>
#include <vector>

namespace Scine {
namespace Utils {
namespace UniversalSettings {

/*!
 * SettingDescriptor for a list of options, of which one must be chosen.
 */
class OptionListDescriptor : public SettingDescriptor {
 public:
  OptionListDescriptor() = default;
  explicit OptionListDescriptor(std::string propertyDescription);

  std::unique_ptr<SettingDescriptor> clone() const override;
  GenericValue getDefaultGenericValue() const override;
  bool validValue(const GenericValue& v) const override;
  std::string explainInvalidValue(const GenericValue& v) const override;

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
  //! Returns the list of option strings
  const std::vector<std::string>& getAllOptions() const;

  /*!
   * @brief Adds a string option to the list of options
   *
   * @param option A new option string
   *
   * @throws OptionAlreadyExistsException if the option already exists
   */
  void addOption(std::string option);
  /**
   * @brief Sets the default option to a particular string
   *
   * @param def The option to select as default
   *
   * @throws OptionDoesNotExistException if the supplied option does not match
   *   an existing option
   */
  void setDefaultOption(const std::string& def);

 private:
  int getIndex(const std::string& option) const;
  std::vector<std::string> options_;
  int defaultIndex_ = 0;
};

inline const std::string& OptionListDescriptor::getDefaultValue() const {
  return getDefaultOption();
}

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */

#endif // UNIVERSALSETTINGS_OPTIONLISTDESCRIPTOR_H
