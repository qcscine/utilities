/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UNIVERSALSETTINGS_GENERICDESCRIPTOR_H
#define UNIVERSALSETTINGS_GENERICDESCRIPTOR_H
/* Internal Headers */
#include "Utils/UniversalSettings/BoolDescriptor.h"
#include "Utils/UniversalSettings/DirectoryDescriptor.h"
#include "Utils/UniversalSettings/DoubleDescriptor.h"
#include "Utils/UniversalSettings/DoubleListDescriptor.h"
#include "Utils/UniversalSettings/FileDescriptor.h"
#include "Utils/UniversalSettings/IntDescriptor.h"
#include "Utils/UniversalSettings/IntListDescriptor.h"
#include "Utils/UniversalSettings/IntListListDescriptor.h"
#include "Utils/UniversalSettings/OptionListDescriptor.h"
#include "Utils/UniversalSettings/StringDescriptor.h"
#include "Utils/UniversalSettings/StringListDescriptor.h"
/* External Headers */
#include <memory>

namespace Scine {
namespace Utils {
namespace UniversalSettings {
/* Forward Declarations */
class DescriptorCollection;
class ParametrizedOptionListDescriptor;
class CollectionListDescriptor;
class GenericValue;

/*!
 * @brief Wrapper around SettingDescriptor that hides the setting type.
 *
 * Makes it possible to handle setting descriptors as objects (but is in principle similar to a pointer to the
 * base class SettingDescriptor).
 */
class GenericDescriptor {
 public:
  //! To store what type the setting descriptor is for
  enum class Type {
    Bool,
    Int,
    Double,
    String,
    File,
    Directory,
    OptionList,
    SettingCollection,
    ParametrizedOptionList,
    IntList,
    IntListList,
    DoubleList,
    StringList,
    CollectionList,
  };

  //!@name Special member functions
  //!@{
  GenericDescriptor() = delete;
  GenericDescriptor(const GenericDescriptor& /*rhs*/);
  GenericDescriptor(GenericDescriptor&&) = default;
  GenericDescriptor& operator=(const GenericDescriptor& /*rhs*/);
  GenericDescriptor& operator=(GenericDescriptor&&) = default;
  ~GenericDescriptor();
  //!@}

  //! Extract what type a setting descriptor is for
  Type getType() const;
  //! Fetches the string describing what the setting descriptor is for
  const std::string& getPropertyDescription() const;
  //! Fetches a type-erased representation of the setting descriptor's default value
  GenericValue getDefaultValue() const;

  //!@name Specific constructors from SettingDescriptor-derived instances
  //!@{
  GenericDescriptor(BoolDescriptor d);
  GenericDescriptor(IntDescriptor d);
  GenericDescriptor(DoubleDescriptor d);
  GenericDescriptor(StringDescriptor d);
  GenericDescriptor(FileDescriptor d);
  GenericDescriptor(DirectoryDescriptor d);
  GenericDescriptor(OptionListDescriptor d);
  GenericDescriptor(DescriptorCollection d);
  GenericDescriptor(ParametrizedOptionListDescriptor d);
  GenericDescriptor(IntListDescriptor d);
  GenericDescriptor(IntListListDescriptor d);
  GenericDescriptor(DoubleListDescriptor d);
  GenericDescriptor(StringListDescriptor d);
  GenericDescriptor(CollectionListDescriptor d);
  //!@}

  //!@name Figure out which SettingsDescriptor derived type is stored
  //!@{
  bool relatesToBool() const;
  bool relatesToInt() const;
  bool relatesToDouble() const;
  bool relatesToString() const;
  bool relatesToFile() const;
  bool relatesToDirectory() const;
  bool relatesToOptionList() const;
  bool relatesToSettingCollection() const;
  bool relatesToParametrizedOptionList() const;
  bool relatesToIntList() const;
  bool relatesToIntListList() const;
  bool relatesToDoubleList() const;
  bool relatesToStringList() const;
  bool relatesToCollectionList() const;
  //!@}

  //!@name Fetch reference to specific descriptor
  //!@{
  const SettingDescriptor& getDescriptor() const;
  const BoolDescriptor& getBoolDescriptor() const;
  const IntDescriptor& getIntDescriptor() const;
  const DoubleDescriptor& getDoubleDescriptor() const;
  const StringDescriptor& getStringDescriptor() const;
  const FileDescriptor& getFileDescriptor() const;
  const DirectoryDescriptor& getDirectoryDescriptor() const;
  const OptionListDescriptor& getOptionListDescriptor() const;
  const DescriptorCollection& getSettingCollectionDescriptor() const;
  const ParametrizedOptionListDescriptor& getParametrizedOptionListDescriptor() const;
  const IntListDescriptor& getIntListDescriptor() const;
  const IntListListDescriptor& getIntListListDescriptor() const;
  const DoubleListDescriptor& getDoubleListDescriptor() const;
  const StringListDescriptor& getStringListDescriptor() const;
  const CollectionListDescriptor& getCollectionListDescriptor() const;
  //!@}

  //!@name General interface
  //!@{
  template<typename T>
  bool is() const {
    // do pointer casting to avoid exceptions
    return dynamic_cast<T*>(descriptor_.get()) != nullptr;
  }

  template<typename T>
  T& get() {
    try {
      return dynamic_cast<T&>(*descriptor_);
    }
    catch (...) {
      throw InvalidDescriptorConversionException(*descriptor_);
    }
  }

  template<typename T>
  const T& get() const {
    try {
      return dynamic_cast<const T&>(*descriptor_);
    }
    catch (...) {
      throw InvalidDescriptorConversionException(*descriptor_);
    }
  }
  //!@}

 private:
  std::unique_ptr<SettingDescriptor> descriptor_;
};

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */

#endif // UNIVERSALSETTINGS_GENERICDESCRIPTOR_H
