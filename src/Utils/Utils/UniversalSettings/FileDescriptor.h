/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UNIVERSALSETTINGS_FILEDESCRIPTOR_H
#define UNIVERSALSETTINGS_FILEDESCRIPTOR_H
/* Internal Headers */
#include "Utils/UniversalSettings/SettingDescriptor.h"
/* External Headers */
#include <string>
#include <vector>

namespace Scine {
namespace Utils {
namespace UniversalSettings {

/*!
 * @brief SettingDescriptor for a file path.
 */
class FileDescriptor : public SettingDescriptor {
 public:
  //! Describes what kind of file the path references
  enum class FileType {
    //! Any kind of file
    Any,
    //! An executable file (permissions set to user-executable)
    Executable
  };

  /**
   * @brief Constructor
   *
   * @param propertyDescription A string describing what the file path is for
   */
  explicit FileDescriptor(std::string propertyDescription);

  std::unique_ptr<SettingDescriptor> clone() const override;
  GenericValue getDefaultGenericValue() const override;
  bool validValue(const GenericValue& v) const override;

  //! Get the default file path
  const std::string& getDefaultValue() const;
  //! Set the default file path
  void setDefaultValue(std::string def);

  //! Get whether the referenced file must already exist
  bool fileMustAlreadyExist() const;
  //! Set whether the referenced file must already exist
  void setFileMustAlreadyExist(bool b);

  //! Get the kind of file referenced
  FileType getFileType() const;
  //! Set what kind of file is referenced
  void setFileType(FileType type);

  /*! Get Qt5 name filter-compatible strings
   *
   * @note Name filters are used for the GUI only and should be compatible with
   * the Qt syntax (see http://doc.qt.io/qt-5/qfiledialog.html#details).
   */
  const std::vector<std::string>& getNameFilters() const;
  /*! Add a Qt5 name filter-compatible string
   *
   * @note Name filters are used for the GUI only and should be compatible with
   * the Qt syntax (see http://doc.qt.io/qt-5/qfiledialog.html#details).
   */
  void addNameFilter(std::string nameFilter);

 private:
  std::string defaultValue_{};
  std::vector<std::string> nameFilters_{};
  bool existingFile_ = true;
  FileType type_ = FileType::Any;
};

inline FileDescriptor::FileDescriptor(std::string propertyDescription)
  : SettingDescriptor(std::move(propertyDescription)) {
}

inline const std::string& FileDescriptor::getDefaultValue() const {
  return defaultValue_;
}

inline void FileDescriptor::setDefaultValue(std::string def) {
  defaultValue_ = std::move(def);
}

inline bool FileDescriptor::fileMustAlreadyExist() const {
  return existingFile_;
}

inline void FileDescriptor::setFileMustAlreadyExist(bool b) {
  existingFile_ = b;
}

inline FileDescriptor::FileType FileDescriptor::getFileType() const {
  return type_;
}

inline void FileDescriptor::setFileType(FileDescriptor::FileType type) {
  type_ = type;
}

inline const std::vector<std::string>& FileDescriptor::getNameFilters() const {
  return nameFilters_;
}

inline void FileDescriptor::addNameFilter(std::string nameFilter) {
  nameFilters_.emplace_back(std::move(nameFilter));
}

inline std::unique_ptr<SettingDescriptor> FileDescriptor::clone() const {
  return std::make_unique<FileDescriptor>(*this);
}

inline GenericValue FileDescriptor::getDefaultGenericValue() const {
  return GenericValue::fromString(getDefaultValue());
}

inline bool FileDescriptor::validValue(const GenericValue& v) const {
  return v.isString();
}

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */

#endif // UNIVERSALSETTINGS_FILEDESCRIPTOR_H
