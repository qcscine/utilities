/**
 * @file FormattedStreamHandler.h
 * @brief Defines interface for classes handling formatted IO from streams
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef CORE_FORMATTED_STREAM_HANDLER_H_
#define CORE_FORMATTED_STREAM_HANDLER_H_

#include <istream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace Scine {
namespace Utils {

// Forward-declarations
class AtomCollection;
class BondOrderCollection;

/**
 * @brief The interface for all classes handling formatted streaming IO
 */
class FormattedStreamHandler {
 public:
  //! Which operations are supported for a particular format
  enum class SupportType { ReadOnly, ReadWrite, WriteOnly };

  //! Pair of a format filetype suffix and a @p SupportType
  using FormatSupportPair = std::pair<std::string, SupportType>;

  //*@name Exception types
  //*@{
  /**
   * @brief Exception thrown if the desired format is unsupported by the
   *   handler
   */
  struct FormatUnsupportedException final : public std::exception {
    const char* what() const noexcept final {
      return "The selected format is not supported.";
    }
  };
  /**
   * @brief Exception thrown if the specified format and the input stream
   *   mismatch
   */
  struct FormatMismatchException final : public std::exception {
    const char* what() const noexcept final {
      return "The selected format and the provided stream do not match.";
    }
  };
  /**
   * @brief Exception thrown if the specified number of atoms was exceeded
   */
  struct AtomNumberMismatchException final : public std::exception {
    const char* what() const noexcept final {
      return "Encountered more or fewer atoms than specified.";
    }
  };
  /**
   * @brief Exception thrown if the function yielding or taking a
   *   BondOrderCollection reads from or writes to a stream whose format does
   *   not include bond order information.
   */
  struct NoBondInformationException final : public std::exception {
    const char* what() const noexcept final {
      return "The stream format does not include bond information.";
    }
  };
  //*@}

  /// @brief Default Destructor.
  virtual ~FormattedStreamHandler() = default;

  /**
   * @brief Reads from an input stream, extracting element types, positions and
   *   bond orders (if present)
   * @param is The input stream to read
   * @throws FormatUnsupportedException If the desired format is unsupported
   * @returns An Atomcollection and BondOrderCollection
   */
  virtual std::pair<Utils::AtomCollection, Utils::BondOrderCollection> read(std::istream& is,
                                                                            const std::string& format) const = 0;

  /**
   * @brief Writes element types and positional information to an arbitrary
   *   filetype
   * @param os The output stream to write to
   * @param format The format to write
   * @param atoms The element type and positional information to write to a file
   * @throws FormatUnsupportedException If the desired format is unsupported
   */
  virtual void write(std::ostream& os, const std::string& format, const Utils::AtomCollection& atoms,
                     const std::string& comment = "") const = 0;

  /**
   * @brief Writes element types and positional information to an arbitrary
   *   filetype
   * @param os The output stream to write to.
   * @param format The format to write
   * @param atoms The element type and positional information to write to a file
   * @param bondOrders The bond order collection to write to file
   * @throws FormatUnsupportedException If the desired filetype is unsupported
   * @throws NoBondInformationException If the specified filetype does not
   *   contain bond order information
   */
  virtual void write(std::ostream& os, const std::string& format, const Utils::AtomCollection& atoms,
                     const Utils::BondOrderCollection& bondOrders, const std::string& comment = "") const = 0;

  /**
   * @brief Returns a list of supported file formats.
   * @note These may, but do not have to be identical to each format's
   *   conventional file suffix
   */
  virtual std::vector<FormatSupportPair> formats() const = 0;

  /**
   * @brief Check whether a particular operation for a format is supported
   * @param format The format filetype suffix to check for
   * @param operation The operation for which support is queried
   * @returns true if the operation is supported for that format
   */
  virtual bool formatSupported(const std::string& format, SupportType operation = SupportType::ReadWrite) const = 0;

  /**
   * @brief Getter for the name of the FormattedStreamHandler
   * @returns Returns the name of the FormattedStreamHandler.
   */
  virtual std::string name() const = 0;
};

} /* namespace Utils */
} /* namespace Scine */

#endif /* CORE_CALCULATOR_H_ */
