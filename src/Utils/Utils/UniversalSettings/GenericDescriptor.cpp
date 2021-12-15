/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
/* Internal Headers */
#include "Utils/UniversalSettings/GenericDescriptor.h"
#include "Utils/UniversalSettings/CollectionListDescriptor.h"
#include "Utils/UniversalSettings/DescriptorCollection.h"
#include "Utils/UniversalSettings/Exceptions.h"
#include "Utils/UniversalSettings/ParametrizedOptionListDescriptor.h"
/* External Headers */
#include <stdexcept>

namespace Scine {
namespace Utils {
namespace UniversalSettings {

namespace {

/* Implement the cast in terms of pointers instead of references to avoid
 * throwing std::bad_cast, which can be annoying while debugging and breaking
 * on all exceptions.
 */
template<typename DerivedClass>
const DerivedClass* convertPtr(const SettingDescriptor* d) {
  const auto* r = dynamic_cast<const DerivedClass*>(d);
  return r;
}

template<typename DerivedClass>
const DerivedClass& convert(const SettingDescriptor& d) {
  const auto* r = convertPtr<DerivedClass>(&d);
  if (r == nullptr) {
    throw InvalidDescriptorConversionException(d);
  }
  return *r;
}

template<typename DerivedClass>
bool verifyConversion(const SettingDescriptor& d) {
  auto ptr = convertPtr<DerivedClass>(&d);
  return ptr != nullptr;
}

} // namespace

GenericDescriptor::~GenericDescriptor() = default;

GenericDescriptor::GenericDescriptor(const GenericDescriptor& rhs) {
  descriptor_ = rhs.descriptor_->clone();
}

GenericDescriptor& GenericDescriptor::operator=(const GenericDescriptor& rhs) {
  descriptor_ = rhs.descriptor_->clone();
  return *this;
}

const std::string& GenericDescriptor::getPropertyDescription() const {
  return descriptor_->getPropertyDescription();
}

GenericValue GenericDescriptor::getDefaultValue() const {
  return descriptor_->getDefaultGenericValue();
}

GenericDescriptor::GenericDescriptor(BoolDescriptor d) {
  descriptor_ = std::make_unique<BoolDescriptor>(std::move(d));
}

GenericDescriptor::GenericDescriptor(IntDescriptor d) {
  descriptor_ = std::make_unique<IntDescriptor>(std::move(d));
}

GenericDescriptor::GenericDescriptor(DoubleDescriptor d) {
  descriptor_ = std::make_unique<DoubleDescriptor>(std::move(d));
}

GenericDescriptor::GenericDescriptor(StringDescriptor d) {
  descriptor_ = std::make_unique<StringDescriptor>(std::move(d));
}

GenericDescriptor::GenericDescriptor(FileDescriptor d) {
  descriptor_ = std::make_unique<FileDescriptor>(std::move(d));
}

GenericDescriptor::GenericDescriptor(DirectoryDescriptor d) {
  descriptor_ = std::make_unique<DirectoryDescriptor>(std::move(d));
}

GenericDescriptor::GenericDescriptor(OptionListDescriptor d) {
  descriptor_ = std::make_unique<OptionListDescriptor>(std::move(d));
}

GenericDescriptor::GenericDescriptor(DescriptorCollection d) {
  descriptor_ = std::make_unique<DescriptorCollection>(std::move(d));
}

GenericDescriptor::GenericDescriptor(ParametrizedOptionListDescriptor d) {
  descriptor_ = std::make_unique<ParametrizedOptionListDescriptor>(std::move(d));
}

GenericDescriptor::GenericDescriptor(IntListDescriptor d) {
  descriptor_ = std::make_unique<IntListDescriptor>(std::move(d));
}

GenericDescriptor::GenericDescriptor(DoubleListDescriptor d) {
  descriptor_ = std::make_unique<DoubleListDescriptor>(std::move(d));
}

GenericDescriptor::GenericDescriptor(StringListDescriptor d) {
  descriptor_ = std::make_unique<StringListDescriptor>(std::move(d));
}

GenericDescriptor::GenericDescriptor(CollectionListDescriptor d) {
  descriptor_ = std::make_unique<CollectionListDescriptor>(std::move(d));
}

GenericDescriptor::Type GenericDescriptor::getType() const {
  if (relatesToBool()) {
    return Type::Bool;
  }
  if (relatesToInt()) {
    return Type::Int;
  }
  if (relatesToDouble()) {
    return Type::Double;
  }
  if (relatesToString()) {
    return Type::String;
  }
  if (relatesToFile()) {
    return Type::File;
  }
  if (relatesToDirectory()) {
    return Type::Directory;
  }
  if (relatesToOptionList()) {
    return Type::OptionList;
  }
  if (relatesToSettingCollection()) {
    return Type::SettingCollection;
  }
  if (relatesToParametrizedOptionList()) {
    return Type::ParametrizedOptionList;
  }
  if (relatesToIntList()) {
    return Type::IntList;
  }
  if (relatesToDoubleList()) {
    return Type::DoubleList;
  }
  if (relatesToStringList()) {
    return Type::StringList;
  }
  if (relatesToCollectionList()) {
    return Type::CollectionList;
  }
  throw Exception("GenericDescriptor has an unknown type.");
}

bool GenericDescriptor::relatesToBool() const {
  return verifyConversion<BoolDescriptor>(*descriptor_);
}

bool GenericDescriptor::relatesToInt() const {
  return verifyConversion<IntDescriptor>(*descriptor_);
}

bool GenericDescriptor::relatesToDouble() const {
  return verifyConversion<DoubleDescriptor>(*descriptor_);
}

bool GenericDescriptor::relatesToString() const {
  return verifyConversion<StringDescriptor>(*descriptor_);
}

bool GenericDescriptor::relatesToFile() const {
  return verifyConversion<FileDescriptor>(*descriptor_);
}

bool GenericDescriptor::relatesToDirectory() const {
  return verifyConversion<DirectoryDescriptor>(*descriptor_);
}

bool GenericDescriptor::relatesToOptionList() const {
  return verifyConversion<OptionListDescriptor>(*descriptor_);
}

bool GenericDescriptor::relatesToSettingCollection() const {
  return verifyConversion<DescriptorCollection>(*descriptor_);
}

bool GenericDescriptor::relatesToParametrizedOptionList() const {
  return verifyConversion<ParametrizedOptionListDescriptor>(*descriptor_);
}

bool GenericDescriptor::relatesToIntList() const {
  return verifyConversion<IntListDescriptor>(*descriptor_);
}

bool GenericDescriptor::relatesToDoubleList() const {
  return verifyConversion<DoubleListDescriptor>(*descriptor_);
}

bool GenericDescriptor::relatesToStringList() const {
  return verifyConversion<StringListDescriptor>(*descriptor_);
}

bool GenericDescriptor::relatesToCollectionList() const {
  return verifyConversion<CollectionListDescriptor>(*descriptor_);
}

const SettingDescriptor& GenericDescriptor::getDescriptor() const {
  return *descriptor_;
}

const BoolDescriptor& GenericDescriptor::getBoolDescriptor() const {
  return convert<BoolDescriptor>(*descriptor_);
}

const IntDescriptor& GenericDescriptor::getIntDescriptor() const {
  return convert<IntDescriptor>(*descriptor_);
}

const DoubleDescriptor& GenericDescriptor::getDoubleDescriptor() const {
  return convert<DoubleDescriptor>(*descriptor_);
}

const StringDescriptor& GenericDescriptor::getStringDescriptor() const {
  return convert<StringDescriptor>(*descriptor_);
}

const FileDescriptor& GenericDescriptor::getFileDescriptor() const {
  return convert<FileDescriptor>(*descriptor_);
}

const DirectoryDescriptor& GenericDescriptor::getDirectoryDescriptor() const {
  return convert<DirectoryDescriptor>(*descriptor_);
}

const OptionListDescriptor& GenericDescriptor::getOptionListDescriptor() const {
  return convert<OptionListDescriptor>(*descriptor_);
}

const DescriptorCollection& GenericDescriptor::getSettingCollectionDescriptor() const {
  return convert<DescriptorCollection>(*descriptor_);
}

const ParametrizedOptionListDescriptor& GenericDescriptor::getParametrizedOptionListDescriptor() const {
  return convert<ParametrizedOptionListDescriptor>(*descriptor_);
}

const IntListDescriptor& GenericDescriptor::getIntListDescriptor() const {
  return convert<IntListDescriptor>(*descriptor_);
}

const DoubleListDescriptor& GenericDescriptor::getDoubleListDescriptor() const {
  return convert<DoubleListDescriptor>(*descriptor_);
}

const StringListDescriptor& GenericDescriptor::getStringListDescriptor() const {
  return convert<StringListDescriptor>(*descriptor_);
}

const CollectionListDescriptor& GenericDescriptor::getCollectionListDescriptor() const {
  return convert<CollectionListDescriptor>(*descriptor_);
}

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */
