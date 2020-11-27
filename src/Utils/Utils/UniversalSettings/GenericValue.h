/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UNIVERSALSETTINGS_GENERICVALUE_H
#define UNIVERSALSETTINGS_GENERICVALUE_H
/* External Headers */
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace Scine {
namespace Utils {
namespace UniversalSettings {
/* Forward Declarations */
class ValueCollection;
struct ParametrizedOptionValue;

/*!
 *
 * @brief Class that uniformly stores multiple types of values
 *
 * Wrapper around some parameter value that hides the actual value.
 * Makes it possible to handle different types as one object.
 */
class GenericValue {
 public:
  //!@name Types
  //!@{
  using IntList = std::vector<int>;
  using StringList = std::vector<std::string>;
  using CollectionList = std::vector<ValueCollection>;
  //!@}

  //!@name Special member functions
  //!@{
  GenericValue(const GenericValue& rhs);
  GenericValue(GenericValue&& rhs) noexcept;
  GenericValue& operator=(const GenericValue& rhs);
  GenericValue& operator=(GenericValue&& rhs) noexcept;
  ~GenericValue();
  //!@}

  //!@name Factory functions
  //!@{
  static GenericValue fromBool(bool v);
  static GenericValue fromInt(int v);
  static GenericValue fromDouble(double v);
  static GenericValue fromString(std::string s);
  static GenericValue fromCollection(ValueCollection s);
  static GenericValue fromOptionWithSettings(ParametrizedOptionValue v);

  static GenericValue fromIntList(IntList v);
  static GenericValue fromStringList(const StringList& v);
  static GenericValue fromCollectionList(CollectionList v);
  //!@}

  //!@name Type determinators
  //!@{
  bool isBool() const;
  bool isInt() const;
  bool isDouble() const;
  bool isString() const;
  bool isCollection() const;
  bool isOptionWithSettings() const;

  bool isIntList() const;
  bool isStringList() const;
  bool isCollectionList() const;
  //!@}

  //!@name Type converters
  //!@{
  bool toBool() const;
  int toInt() const;
  double toDouble() const;
  std::string toString() const;
  ValueCollection toCollection() const;
  ParametrizedOptionValue toOptionWithSettings() const;

  IntList toIntList() const;
  StringList toStringList() const;
  CollectionList toCollectionList() const;
  //!@}

 private:
  /*!
   * @brief Private constructor
   *
   * This is to avoid instantiation with an empty boost::any member where no
   * conversions are possible and all type determinators fail.
   */
  GenericValue();
  struct Impl; // Use pimpl idiom to hide boost::any dependency.
  using ImplPtr = std::unique_ptr<Impl>;
  ImplPtr pimpl_;
};

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */

#endif // UNIVERSALSETTINGS_GENERICVALUE_H
