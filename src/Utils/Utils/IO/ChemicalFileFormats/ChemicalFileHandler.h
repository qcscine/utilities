/**
 * @file
 * @brief File-level I/O on chemical file formats
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef INCLUDE_SCINE_UTILS_FILE_HANDLER_H
#define INCLUDE_SCINE_UTILS_FILE_HANDLER_H

#include <Utils/Bonds/BondOrderCollection.h>
#include <stdexcept>
#include <string>

namespace Scine {
namespace Utils {

// Forward-declarations
class AtomCollection;

/**
 * @brief Handles the reading and writing of chemical file formats.
 *
 * Matches suffixes to chemical file formats and then dispatches the calls to
 * appropriate streaming I/O handlers. Catches I/O errors.
 */
class ChemicalFileHandler {
 public:
  /**
   * @brief Exception thrown if the file does not exist or cannot be opened
   */
  struct FileInaccessibleException final : public std::exception {
    const char* what() const noexcept final {
      return "The specified file is inaccessible";
    }
  };

  /**
   * @brief Reads an atom collection and tries to read a bond order collection
   *   from a file
   *
   * @param filename The path to the file to read
   *
   * Supported file formats are currently:
   * - mol
   * - xyz
   *
   * If `obabel` is found in your path, all formats that openbabel supports are
   * also available.
   *
   * @note The format of a file is deduced from its suffix.
   * @throws FileInaccessibleException if the file does not exist or cannot be
   *   opened
   * @throws FormatUnsupportedException if the file suffix cannot be matched
   *   to a supported file format
   * @throws FormatMismatchException if the file's format does not match the
   *   parser's expectations (e.g. if the suffix does not match the format)
   * @throws NoBondInformationException if the format does not contain bond
   *   order information
   *
   * @returns an AtomCollection and a BondOrderCollection read from the file.
   *   The BondOrderCollection may be empty if the file format does not include
   *   such information
   */
  static std::pair<AtomCollection, BondOrderCollection> read(const std::string& filename);

  /**
   * @brief Writes an atom collection to a file
   *
   * Supported file formats are currently
   * - xyz
   * - mol
   *
   * @param filename The path to the file to write
   * @param atoms The atom collection to commit to a file
   * @param comment An optional comment to insert into the written file
   *
   * @throws FormatUnsupportedException if the file suffix cannot be matched
   *   to a supported file format
   * @throws FileInaccessibleException if the file cannot be created
   *
   * @note Overwrites existing files
   */
  static void write(const std::string& filename, const AtomCollection& atoms, const std::string& comment = "");

  /**
   * @brief Writes an atom collection and a bond order collection to a file
   *
   * Supported file formats are currently
   * - mol
   *
   * @param filename The path to the file to write
   * @param atoms The atom collection to commit to a file
   * @param comment An optional comment to insert into the written file
   *
   * @throws FormatUnsupportedException if the file suffix cannot be matched
   *   to a supported file format
   * @throws FileInaccessibleException if the file cannot be created
   *
   * @note Overwrites existing files
   */
  static void write(const std::string& filename, const AtomCollection& atoms, const BondOrderCollection& bondOrders,
                    const std::string& comment = "");
};

} // namespace Utils
} // namespace Scine

#endif
