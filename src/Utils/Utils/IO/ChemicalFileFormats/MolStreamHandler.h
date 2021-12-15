/**
 * @file
 * @brief Mol-formatted streaming IO
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef INCLUDE_SCINE_UTILS_IO_MOL_STREAM_HANDLER_H
#define INCLUDE_SCINE_UTILS_IO_MOL_STREAM_HANDLER_H

#include "Utils/IO/ChemicalFileFormats/FormattedStreamHandler.h"
#include "boost/optional/optional_fwd.hpp"

namespace Scine {
namespace Utils {

class AtomCollection;
class BondOrderCollection;

/**
 * @brief Handler for MOL stream IO
 *
 * This class implements the FormattedStreamHandler interface for use with
 * Core's module management. The core stream IO functions are available as
 * static functions, too.
 *
 * @todo Implement V3000 reading/writing
 */
class MolStreamHandler : public FormattedStreamHandler {
 public:
  //! Model string identifier for use in the Module system
  static constexpr const char* model = "MolStreamHandler";

  //!@name Static members
  //!@{
  /**
   * @brief Writes to an output stream in MOL format
   *
   * @param os The stream to write to
   * @param atoms The element and positional data to write
   * @param bondOrdersOption An optional BondOrderCollection
   * @param formatVersion The string (matching an entry in versionStrings)
   *   specifying which MOL format version should be written
   * @param comment An optional comment to insert into the written file
   *
   * @note Currently only implements the V2000 format.
   */
  static void write(std::ostream& os, const AtomCollection& atoms, const boost::optional<BondOrderCollection>& bondOrdersOption,
                    const std::string& formatVersion, const std::string& comment = "");

  /**
   * @brief Reads from an input stream in MOL format
   *
   * @param is The stream to read from
   *
   * @return Elemental and positional data. Optionally also bond order
   *   information.
   */
  static std::pair<AtomCollection, BondOrderCollection> read(std::istream& is);
  //!@}

  //!@name FormattedStreamHandler interface
  //!@{
  std::pair<AtomCollection, BondOrderCollection> read(std::istream& is, const std::string& format) const final;

  void write(std::ostream& os, const std::string& format, const AtomCollection& atoms, const std::string& comment = "") const final;

  void write(std::ostream& os, const std::string& format, const AtomCollection& atoms,
             const BondOrderCollection& bondOrders, const std::string& comment = "") const final;

  std::vector<FormatSupportPair> formats() const final;

  bool formatSupported(const std::string& format, SupportType /* operation */) const final;

  std::string name() const final;
  //!@}
};

} // namespace Utils
} // namespace Scine

#endif
