/**
 * @file
 * @brief Pdb-formatted streaming IO
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef INCLUDE_SCINE_UTILS_IO_PDB_STREAM_HANDLER_H
#define INCLUDE_SCINE_UTILS_IO_PDB_STREAM_HANDLER_H

#include "Utils/IO/ChemicalFileFormats/FormattedStreamHandler.h"
#include "boost/optional/optional_fwd.hpp"

namespace Scine {
namespace Utils {

struct PdbFileData {
  // The content of the PDB file
  std::string content;
  // The header of the PDB file including diverse information.
  std::string header;
  // The atom block of the PDB file
  std::string atomBlock;
  // The bond order block of the PDB file
  std::string connectivityBlock;
  // A vector of overlaying substructure identifiers
  std::vector<std::string> overlayIdentifiers;
  // The number of Atoms in the file
  int nAtoms = 0;
};

class AtomCollection;
class BondOrderCollection;

/**
 * @class PdbStreamHandler PdbStreamHandler.h
 * @brief A class that handles the reading and writing of PDB file formats
 *
 * @note Since the connectivity is only represented in some PDB's and usually very error-prone, it
 * is not parsed.
 */
class PdbStreamHandler : public FormattedStreamHandler {
 public:
  //! Model string identifier for use in the Module system
  static constexpr const char* model = "PdbStreamHandler";
  //!@name Static members
  //!@{
  /**
   * @brief Writes to an output stream in PDB format
   *
   * @param os The stream to write to
   * @param atoms The element and positional data to write
   *
   * @note Currently only writes "UNX" as residue identifier.
   */
  static void write(std::ostream& os, const AtomCollection& atoms, const std::string& comment = "");
  /**
   * @brief Reads from an input stream in PDB format
   *
   * @param is The stream to read from
   *
   * @return Elemental and positional data.
   * @note Reading of hydrogens or solvent can be set to false via setReadH(false) or setReadHOH(false)
   */
  std::vector<AtomCollection> read(std::istream& is) const;
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
  //!@
  // Function to set whether hydrogen atoms should be parse
  void setReadH(bool includeH);
  // Function to set whether solvent molecules should be parsed
  void setReadHOH(bool includeHOH);
  void setSubstructureID(int substructureID);

 private:
  static void extractContent(std::istream& is, PdbFileData& data);
  bool includeH_ = true;
  bool includeHOH_ = true;
  unsigned int substructureID_ = 0;
};

} // namespace Utils
} // namespace Scine

#endif // INCLUDE_SCINE_UTILS_IO_PDB_STREAM_HANDLER_H
