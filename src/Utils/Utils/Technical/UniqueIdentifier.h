/**
 * @file
 * @brief A file containing definitions of classes that are just different names
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_UNIQUEIDENTIFIER_H
#define UTILS_UNIQUEIDENTIFIER_H

#include <memory>
#include <string>

namespace Scine {
namespace Utils {

/**
 * @class UniqueIdentifier UniqueIdentifier.h
 * @brief   Class for a unique identifier (handle).
 *          It can f.i. be used to identify an instance unequivocally.
 *          Uses pimpl idiom to hide boost::uuid dependency.
 */
class UniqueIdentifier {
 public:
  /** @brief Default Constructor */
  UniqueIdentifier();
  /** @brief Default Destructor */
  ~UniqueIdentifier();
  /** @brief Constructor that takes in a unique ID as a const reference */
  UniqueIdentifier(const UniqueIdentifier& /*rhs*/);
  /** @brief Constructor that takes in a unique ID as an rvalue reference */
  UniqueIdentifier(UniqueIdentifier&& /*rhs*/) noexcept;
  /** @brief Assignment operator */
  UniqueIdentifier& operator=(const UniqueIdentifier& /*rhs*/);
  /** @brief Assignment operator */
  UniqueIdentifier& operator=(UniqueIdentifier&& /*rhs*/) noexcept;

  /**
   * @brief Return a string representation of the unique identifier.
   */
  std::string getStringRepresentation() const;

  /**
   * @brief Implementation of the "equal" operator.
   */
  bool operator==(const UniqueIdentifier& rhs) const;
  /**
   * @brief Implementation of the "not equal" operator.
   */
  bool operator!=(const UniqueIdentifier& rhs) const;
  /**
   * @brief Implementation of the "smaller than" operator.
   */
  bool operator<(const UniqueIdentifier& rhs) const; // to allow use as key in std::map

 private:
  /* pimpl idiom */
  struct Impl; // pimpl idiom
  std::unique_ptr<Impl> pimpl;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_UNIQUEIDENTIFIER_H
