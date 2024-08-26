/**
 * @file
 * @brief Pdb-formatted streaming IO
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef INCLUDE_SCINE_UTILS_IO_PDB_STREAM_HANDLER_H
#define INCLUDE_SCINE_UTILS_IO_PDB_STREAM_HANDLER_H

#include "Utils/Bonds/BondOrderCollection.h"
#include "Utils/IO/ChemicalFileFormats/FormattedStreamHandler.h"
#include "boost/optional/optional_fwd.hpp"

namespace Scine {
namespace Utils {

class AtomCollection;
class MolecularTrajectory;
class BondOrderCollection;
struct PdbFileData {
  std::vector<AtomCollection> atomCollections;
  BondOrderCollection bondOrderCollection;
};

class Atom;
using ResidueInformation = std::tuple<std::string, std::string, std::string, unsigned int>;

/**
 * @class PdbStreamHandler PdbStreamHandler.h
 * @brief A class that handles the reading and writing of PDB file formats
 *
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
   * @param bondOrders Bond orders to write to the file.
   * @param comment Optional comment for the file.
   * @param trajectoryFormat If true, statements that end the pdb file are omitted.
   * @param counter Counter for the model entry.
   *
   * @note Currently only writes "UNX" as residue identifier.
   */
  static void write(std::ostream& os, const AtomCollection& atoms,
                    const BondOrderCollection& bondOrders = BondOrderCollection(), const std::string& comment = "",
                    bool trajectoryFormat = false, unsigned int counter = 0);
  /**
   * @brief Write a pdb trajectory.
   * @param os         The output stream.
   * @param m          The molecular trajectory.
   * @param bondOrders Optional bond orders.
   * @param comment    Optional comment.
   */
  static void writeTrajectory(std::ostream& os, const MolecularTrajectory& m,
                              const BondOrderCollection& bondOrders = BondOrderCollection(), const std::string& comment = "");
  /**
   * @brief Reads from an input stream in PDB format
   *
   * @param is The stream to read from
   *
   * @return Elemental and positional data.
   * @note Reading of hydrogen atoms can be set to false via setReadH(false)
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
  // Function to set whether hydrogen atoms should be parsed
  void setReadH(bool includeH);
  // Function to set whether solvent molecules should be parsed
  void setSubstructureID(int substructureID);
  ///@brief Convert the string encoding the residue sequence number to an int.
  static unsigned int sequenceNumberStringToInt(const std::string& sequenceNumberString);
  static std::string sequenceNumberIntToString(const unsigned int& sequenceInteger);

 private:
  ///@brief Create an atom from the coordinate information in a PDB file line.
  inline static Atom getAtomFromPdbLine(const std::string& line);
  ///@brief Get the atom-like index for a terminus definition in a PDB file line.
  inline static unsigned int getTerminusIndexFromLine(const std::string& line);
  ///@brief Create a ResidueInformation object from an ATOM/HETATM PDB file line.
  inline static ResidueInformation getResidueInformationFromPdbLine(const std::string& line);
  ///@brief Remove all spaces from the given string.
  inline static std::string removeAllSpacesFromString(std::string string);
  /**
   * @brief To make our life more difficult, termini can be encoded in the pdb file. These termini count
   * towards the index used in the CONECT statements. Therefore, we must adjust the atom indices
   * in the CONECT statements by the number of termini with an index lower than this index. This function
   * performs this shift.
   * @param triplets The bonding information encoded as triplets, i.e., the bond order is encoded as a triplet of atom
   *                 index, atom index, bond order.
   * @param termini The atom-like indices of the terminis defined in the PDB file.
   */
  inline static void shiftBondIndicesByTermini(std::vector<Eigen::Triplet<double>>& triplets,
                                               const std::vector<unsigned int>& termini);
  ///@brief Shift the atom index by the number of termini with lower atom-like index (see shiftBondIndicesByTermini for
  /// more information).
  inline static unsigned int shiftAtomIndexByTermini(unsigned int iAtom, const std::vector<unsigned int>& termini);
  ///@brief Parse the PDB file.
  PdbFileData extractContent(std::istream& is, bool onlyOneStructure = false) const;

  ///@brief True if the line encodes ATOM/HETATM information
  static bool isAtomLine(const std::string& line);
  ///@brief True if the line signals the end of a molecular model.
  static bool isEndModelLine(const std::string& line);
  ///@brief True if the line is part of the CONECT information.
  static bool isConnectLine(const std::string& line);
  ///@brief True if the line encodes a terminus.
  static bool isTerminusLine(const std::string& line);
  bool includeH_ = true;
  unsigned int substructureID_ = 0;
  std::vector<unsigned int> residueIndices_;
  ///@brief Extract bonding information from a CONECT line.
  static std::vector<Eigen::Triplet<double>> getConnectTriplet(const std::string& line);
};

} // namespace Utils
} // namespace Scine

#endif // INCLUDE_SCINE_UTILS_IO_PDB_STREAM_HANDLER_H
