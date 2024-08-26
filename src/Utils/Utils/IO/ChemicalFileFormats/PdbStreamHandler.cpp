/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/IO/ChemicalFileFormats/PdbStreamHandler.h"
#include "Utils/Constants.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Geometry/ElementInfo.h"
#include "Utils/MolecularTrajectory.h"
#include <fstream>
#include <iomanip>
#include <iostream>

namespace Scine {
namespace Utils {

std::pair<AtomCollection, BondOrderCollection> PdbStreamHandler::read(std::istream& is, const std::string& format) const {
  if (format != "pdb") {
    throw FormattedStreamHandler::FormatUnsupportedException();
  }
  PdbFileData data = extractContent(is, true);
  if (data.atomCollections.size() < substructureID_) {
    throw std::runtime_error("Structure index out of range.");
  }
  return std::make_pair(data.atomCollections[substructureID_], data.bondOrderCollection);
}

void PdbStreamHandler::write(std::ostream& os, const std::string& format, const AtomCollection& atoms,
                             const std::string& comment) const {
  if (format != "pdb") {
    throw FormattedStreamHandler::FormatUnsupportedException();
  }
  write(os, atoms, BondOrderCollection(), comment);
}

void PdbStreamHandler::write(std::ostream& os, const std::string& format, const AtomCollection& atoms,
                             const BondOrderCollection& bondOrders, const std::string& comment) const {
  if (format != "pdb") {
    throw FormattedStreamHandler::FormatUnsupportedException();
  }
  write(os, atoms, bondOrders, comment);
}

std::vector<PdbStreamHandler::FormatSupportPair> PdbStreamHandler::formats() const {
  return {{{"pdb", SupportType::ReadWrite}}};
}

bool PdbStreamHandler::formatSupported(const std::string& format, SupportType /* operation */
) const {
  return format == "pdb";
}

std::string PdbStreamHandler::name() const {
  return PdbStreamHandler::model;
}

std::vector<Eigen::Triplet<double>> PdbStreamHandler::getConnectTriplet(const std::string& line) {
  std::vector<Eigen::Triplet<double>> triplets;
  std::string connect;
  std::string iString;
  std::string jString;
  std::stringstream lineStream(line);
  lineStream >> connect;
  lineStream >> iString;
  const unsigned int iAtom = std::stoi(iString) - 1;
  while (lineStream >> jString) {
    const unsigned int jAtom = std::stoi(jString) - 1;
    triplets.emplace_back(iAtom, jAtom, 1.0);
  }
  return triplets;
}

unsigned int PdbStreamHandler::getTerminusIndexFromLine(const std::string& line) {
  std::stringstream lineStream(line);
  std::string ter;
  unsigned int index = 0;
  lineStream >> ter >> index;
  return index - 1;
}

PdbFileData PdbStreamHandler::extractContent(std::istream& is, bool onlyOneStructure) const {
  PdbFileData data;
  AtomCollection nextAtomCollection;
  ResidueCollection nextResidueCollection;
  std::string line;
  std::vector<Eigen::Triplet<double>> triplets;
  std::vector<unsigned int> termini;
  bool skipAtoms = false;
  bool skipTermini = false;
  while (std::getline(is, line)) {
    if (isConnectLine(line)) {
      const auto lineTriplets = getConnectTriplet(line);
      triplets.insert(triplets.end(), lineTriplets.begin(), lineTriplets.end());
    }
    else if (!skipAtoms && isAtomLine(line)) {
      const Atom atom = getAtomFromPdbLine(line);
      if (!this->includeH_ && atom.getElementType() == ElementType::H) {
        continue;
      }
      nextResidueCollection.push_back(getResidueInformationFromPdbLine(line));
      nextAtomCollection.push_back(atom);
    }
    else if (!skipTermini && isTerminusLine(line)) {
      termini.push_back(getTerminusIndexFromLine(line));
    }
    else if (isEndModelLine(line)) {
      nextAtomCollection.setResidues(nextResidueCollection);
      data.atomCollections.push_back(nextAtomCollection);
      nextAtomCollection = AtomCollection();
      nextResidueCollection = ResidueCollection();
      skipTermini = true;
      if (onlyOneStructure && data.atomCollections.size() > this->substructureID_) {
        skipAtoms = true;
      }
    }
  }
  if (nextAtomCollection.size() != 0) {
    nextAtomCollection.setResidues(nextResidueCollection);
    data.atomCollections.push_back(nextAtomCollection);
  }
  if (data.atomCollections.empty()) {
    throw std::runtime_error("No atom coordinates found in pdb file.");
  }
  // We cannot parse the CONECT block if we skipped atoms before.
  if (this->includeH_) {
    const unsigned int nAtoms = data.atomCollections[0].size();
    /*
     * To make our life more difficult, termini can be encoded in the pdb file. These termini count
     * towards the index used in the CONECT statements. Therefore, we must adjust the atom indices
     * in the CONECT statements by the number of termini with an index lower than this index.
     */
    shiftBondIndicesByTermini(triplets, termini);
    Eigen::SparseMatrix<double> bondOrderMatrix(nAtoms, nAtoms);
    bondOrderMatrix.setFromTriplets(triplets.begin(), triplets.end());
    BondOrderCollection bondOrderCollection(nAtoms);
    bondOrderCollection.setMatrix(bondOrderMatrix);
    data.bondOrderCollection = bondOrderCollection;
  }
  return data;
}

inline std::string PdbStreamHandler::removeAllSpacesFromString(std::string string) {
  string.erase(std::remove(string.begin(), string.end(), ' '), string.end());
  return string;
}

unsigned int PdbStreamHandler::sequenceNumberStringToInt(const std::string& sequenceNumberString) {
  /*
   * Sometimes the sequence number is given as a hexadecimal number. However, the first 9999 indices
   * are usually still given in decimal. Only then hexadecimal numbers are used. In that case, a
   * hexadecimal number will always contain at least one non-digit character, e.g., A000, A001.
   */
  bool isDecimal = true;
  for (const auto& c : sequenceNumberString) {
    if (!std::isdigit(c)) {
      isDecimal = false;
      break;
    }
  }
  if (isDecimal) {
    return std::stoi(sequenceNumberString);
  }
  else {
    unsigned int x = 0;
    std::stringstream stringstream;
    stringstream << std::hex << sequenceNumberString;
    stringstream >> x;
    return x;
  }
}

inline ResidueInformation PdbStreamHandler::getResidueInformationFromPdbLine(const std::string& line) {
  try {
    const std::string residueNameStr = removeAllSpacesFromString(line.substr(17, 3));
    const std::string atomTypeStr = removeAllSpacesFromString(line.substr(12, 4));
    const std::string chainIdentifier = removeAllSpacesFromString(line.substr(21, 1));
    const std::string residueSequenceNumberStr = removeAllSpacesFromString(line.substr(22, 4));
    unsigned int residueSequenceNumber = 1;
    if (!residueSequenceNumberStr.empty()) {
      residueSequenceNumber = sequenceNumberStringToInt(residueSequenceNumberStr);
    }
    return {residueNameStr, atomTypeStr, chainIdentifier, residueSequenceNumber};
  }
  catch (...) {
    throw std::runtime_error("Unable to read residue information from pdb file.\n"
                             "The problematic line is:\n" +
                             line);
  }
}

inline Atom PdbStreamHandler::getAtomFromPdbLine(const std::string& line) {
  std::string elementStr = removeAllSpacesFromString(line.substr(76, 3));
  elementStr.erase(remove_if(elementStr.begin(), elementStr.end(), [](char c) { return !isalpha(c); }), elementStr.end());
  // Make sure capitalization matches our variant
  std::transform(std::begin(elementStr), std::begin(elementStr) + 1, std::begin(elementStr), ::toupper);
  // Make other letters lowercase
  std::transform(std::begin(elementStr) + 1, std::end(elementStr), std::begin(elementStr) + 1, ::tolower);
  try {
    const ElementType elementType = ElementInfo::elementTypeForSymbol(elementStr);
    const double xCoord = std::stod(line.substr(31, 8));
    const double yCoord = std::stod(line.substr(39, 8));
    const double zCoord = std::stod(line.substr(47, 8));
    Position position(xCoord, yCoord, zCoord);
    position *= Constants::bohr_per_angstrom;
    return Atom(elementType, position);
  }
  catch (...) {
    throw std::runtime_error("Unable to read atom information from pdb file.\n"
                             "The problematic line is:\n" +
                             line);
  }
}

std::vector<AtomCollection> PdbStreamHandler::read(std::istream& is) const {
  return extractContent(is).atomCollections;
}

void PdbStreamHandler::write(std::ostream& os, const AtomCollection& atoms, const BondOrderCollection& bondOrders,
                             const std::string& comment, bool trajectoryFormat, unsigned int counter) {
  if (atoms.getPositions().array().abs().maxCoeff() * Constants::angstrom_per_bohr > 1e+5) {
    throw std::runtime_error("PDB files cannot encode structures with coordinates larger than 1e+5 accurately.\n"
                             "Please use a different file format or ensure that the absolute coordinate values\n"
                             "are small.");
  }
  const unsigned int N = atoms.size();
  const auto& atomResidues = atoms.getResidues();
  if (!comment.empty()) {
    os << "REMARK " << comment << "\n";
  }
  os << "MODEL        " << counter << "\n";
  const PositionCollection positions = atoms.getPositions() * Constants::angstrom_per_bohr;
  for (unsigned int i = 0; i < N; ++i) {
    const std::string element = ElementInfo::symbol(atoms.getElement(int(i)));
    Eigen::Vector3d position = positions.row(i);
    const auto& resLabel = std::get<0>(atomResidues[i]);
    const auto& atomType = std::get<1>(atomResidues[i]);
    const auto& chainLabel = std::get<2>(atomResidues[i]);
    const unsigned& resIndex = std::get<3>(atomResidues[i]);
    const std::string resIndexString = sequenceNumberIntToString(resIndex);
    if (chainLabel.size() > 1) {
      throw std::runtime_error("Chain labels in pdb files may only be one character long.");
    }
    os << "ATOM  " << std::setw(5) << std::right << i + 1 << " " << std::setw(4) << std::left << atomType << " "
       << std::setw(3) << std::right << resLabel << " " << std::setw(1) << std::right << chainLabel << resIndexString
       << "    " << std::setw(8) << std::right << std::fixed << std::setprecision(3) << position(0) << std::setw(8)
       << std::right << std::fixed << std::setprecision(3) << position(1) << std::setw(8) << std::right << std::fixed
       << std::setprecision(3) << position(2) << std::setw(6) << std::right << std::fixed << std::setprecision(2) << 1.0
       << std::setw(6) << std::right << std::fixed << std::setprecision(2) << 0.0 << std::setw(12) << std::right
       << element << "\n";
  }
  os << "ENDMDL" << std::endl;
  if (!trajectoryFormat) {
    if (!bondOrders.empty()) {
      for (unsigned int i = 0; i < N; ++i) {
        const auto bondPartners = bondOrders.getBondPartners(int(i));
        if (!bondPartners.empty()) {
          os << "CONECT" << std::setw(5) << std::right << i + 1;
          for (const auto& j : bondPartners) {
            os << std::setw(5) << std::right << j + 1;
          }
          os << "\n";
        }
      }
    }
    os << "MASTER        0";
    for (unsigned i = 0; i < 7; ++i) {
      os << std::setw(5) << std::right << 0;
    }
    os << std::setw(5) << std::right << N;
    os << std::setw(5) << std::right << 0;
    os << std::setw(5) << std::right << N;
    os << std::setw(5) << std::right << 0 << "\n";
    os << "END" << std::endl;
  }
}

void PdbStreamHandler::setReadH(bool includeH) {
  includeH_ = includeH;
}

void PdbStreamHandler::setSubstructureID(int substructureID) {
  substructureID_ = substructureID;
}

bool PdbStreamHandler::isAtomLine(const std::string& line) {
  return (line.rfind("ATOM", 0) == 0) || (line.rfind("HETATM", 0) == 0);
}

bool PdbStreamHandler::isConnectLine(const std::string& line) {
  return line.find("CONECT") != std::string::npos;
}

bool PdbStreamHandler::isEndModelLine(const std::string& line) {
  return line.find("ENDMDL") != std::string::npos;
  ;
}

void PdbStreamHandler::shiftBondIndicesByTermini(std::vector<Eigen::Triplet<double>>& triplets,
                                                 const std::vector<unsigned int>& termini) {
  if (termini.empty()) {
    return;
  }
  for (auto& triplet : triplets) {
    const unsigned int iAtom = triplet.row();
    const unsigned int jAtom = triplet.col();
    triplet = Eigen::Triplet<double>(shiftAtomIndexByTermini(iAtom, termini), shiftAtomIndexByTermini(jAtom, termini), 1.0);
  }
}

unsigned int PdbStreamHandler::shiftAtomIndexByTermini(unsigned int iAtom, const std::vector<unsigned int>& termini) {
  const unsigned int lessThanIAtom = std::lower_bound(termini.begin(), termini.end(), iAtom) - termini.begin();
  if (iAtom < lessThanIAtom) {
    throw std::runtime_error("Terminus indices in pdb file have incorrect numbering");
  }
  return iAtom - lessThanIAtom;
}

bool PdbStreamHandler::isTerminusLine(const std::string& line) {
  return line.find("TER") != std::string::npos;
  ;
}
void PdbStreamHandler::writeTrajectory(std::ostream& os, const MolecularTrajectory& m,
                                       const BondOrderCollection& bondOrders, const std::string& comment) {
  if (m.empty())
    return;
  if (m.getResidues().empty()) {
    throw std::runtime_error("The residue information is required to write a pdb trajectory file");
  }
  for (int i = 0; i < m.size() - 1; ++i) {
    AtomCollection atoms(m.getElementTypes(), m.at(i));
    atoms.setResidues(m.getResidues());
    write(os, atoms, bondOrders, comment, true, i);
  }
  AtomCollection atoms(m.getElementTypes(), m.at(m.size() - 1));
  atoms.setResidues(m.getResidues());
  write(os, atoms, bondOrders, comment, m.size() - 1);
}
std::string PdbStreamHandler::sequenceNumberIntToString(const unsigned int& sequenceInteger) {
  std::stringstream stream;
  stream << std::setw(4) << std::right;
  if (sequenceInteger > 9999) {
    stream << std::hex;
  }
  if (sequenceInteger > sequenceNumberStringToInt("eeee")) {
    throw std::runtime_error(
        "Residue sequence number is larger than eeee. This cannot be handled while writing pdb files.");
  }
  stream << sequenceInteger;
  std::string result = stream.str();
  result.erase(std::remove_if(result.begin(), result.end(), [](char c) { return c == "'"[0] || c == ','; }), result.end());
  return result;
}

} // namespace Utils
} // namespace Scine
