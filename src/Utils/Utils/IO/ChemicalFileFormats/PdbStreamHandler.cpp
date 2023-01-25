/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/IO/ChemicalFileFormats/PdbStreamHandler.h"
#include "Utils/Bonds/BondOrderCollection.h"
#include "Utils/Constants.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Geometry/ElementInfo.h"
#include <fstream>
#include <iomanip>
#include <iostream>

namespace Scine {
namespace Utils {

std::pair<AtomCollection, BondOrderCollection> PdbStreamHandler::read(std::istream& is, const std::string& format) const {
  if (format != "pdb") {
    throw FormattedStreamHandler::FormatUnsupportedException();
  }
  auto data = read(is);

  if (substructureID_ > data.size()) {
    std::string message = "Cannot parse substructure " + std::to_string(substructureID_) +
                          "when structure size is: " + std::to_string(data.size());
    throw std::runtime_error(message);
  }

  return std::make_pair(data[substructureID_], BondOrderCollection());
}

void PdbStreamHandler::write(std::ostream& os, const std::string& format, const AtomCollection& atoms,
                             const std::string& comment) const {
  if (format != "pdb") {
    throw FormattedStreamHandler::FormatUnsupportedException();
  }

  write(os, atoms, comment);
}

void PdbStreamHandler::write(std::ostream& /* os */, const std::string& format, const AtomCollection& /* atoms */,
                             const BondOrderCollection& /* bondOrders */, const std::string& /* comment */) const {
  if (format != "pdb") {
    throw FormattedStreamHandler::FormatUnsupportedException();
  }

  throw FormattedStreamHandler::NoBondInformationException();
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

// Extract Content from the file
void PdbStreamHandler::extractContent(std::istream& is, PdbFileData& data) {
  is.exceptions(std::ifstream::failbit | std::ifstream::badbit);

  data.content = std::string(std::istreambuf_iterator<char>{is}, {});

  std::istringstream in(data.content);
  std::string line;
  int numModels = 0;
  bool structureHasModels = false;
  while (std::getline(in, line)) {
    // Extract header
    if (!isAtomLine(line) && !isModelLine(line) && !(line.rfind("ENDMDL", 0) == 0)) {
      std::istringstream headerLine(line);
      data.header.append(std::string(std::istreambuf_iterator<char>{headerLine}, {}));
      data.header.append("\n");
    }

    else if (isModelLine(line)) {
      structureHasModels = true;
      numModels++;
    }
    // Extract connectivity block
    else if (line.rfind("CONECT", 0) == 0) {
      std::istringstream connectivityLine(line);
      data.connectivityBlock.append(std::string(std::istreambuf_iterator<char>{connectivityLine}, {}));
      data.connectivityBlock.append("\n");
    }
  }
  if (numModels > 1) {
    data.numModels = numModels;
  }

  in.clear();
  in.seekg(0);

  // iterate over the file again to extract the atom blocks
  if (structureHasModels) {
    data.atomBlocks.resize(data.numModels);
    extractModels(in, line, data);
  }
  else {
    data.atomBlocks.resize(1);
    data.nAtoms = 0;
    extractStructure(in, line, data);
  }
}

std::vector<AtomCollection> PdbStreamHandler::read(std::istream& is) const {
  PdbFileData data;
  extractContent(is, data);

  std::vector<AtomCollection> structures;
  // Create a single identifier if none is parsed
  if (data.overlayIdentifiers.empty()) {
    data.overlayIdentifiers.emplace_back("A");
  }

  auto removeAllSpacesFromString = [](std::string a) -> std::string {
    a.erase(std::remove(a.begin(), a.end(), ' '), a.end());
    return a;
  };

  std::string line;
  for (const auto& atomBlock : data.atomBlocks) {
    for (const auto& s : data.overlayIdentifiers) {
      std::istringstream in(atomBlock);
      AtomCollection structure;
      for (int i = 0; i < data.nAtoms; ++i) {
        std::getline(in, line);
        std::istringstream iss(line);

        if (iss.str().empty()) {
          continue;
        }

        // Get the elements
        std::string elementStr = removeAllSpacesFromString(iss.str().substr(76, 3));
        // get the residue names
        std::string residueNameStr = removeAllSpacesFromString(iss.str().substr(17, 3));
        std::string overlayIdentifier = removeAllSpacesFromString(iss.str().substr(16, 1));

        if (elementStr == "H" && !includeH_) {
          continue;
        }

        if (elementStr == "HOH" && !parseOnlySolvent_) {
          continue;
        }

        elementStr.erase(remove_if(elementStr.begin(), elementStr.end(), [](char c) { return !isalpha(c); }),
                         elementStr.end());
        // Make sure capitalization matches our variant
        std::transform(std::begin(elementStr), std::begin(elementStr) + 1, std::begin(elementStr), ::toupper);
        // Make other letters lowercase
        std::transform(std::begin(elementStr) + 1, std::end(elementStr), std::begin(elementStr) + 1, ::tolower);
        ElementType f;

        try {
          f = ElementInfo::elementTypeForSymbol(elementStr);
        }
        catch (...) {
          throw FormattedStreamHandler::FormatMismatchException();
        }

        // Get the positions
        double x, y, z;
        try {
          x = std::stod(iss.str().substr(31, 8));
          y = std::stod(iss.str().substr(39, 8));
          z = std::stod(iss.str().substr(47, 8));
        }
        catch (...) {
          throw FormattedStreamHandler::FormatMismatchException();
        }

        Position position(x, y, z);
        position *= Constants::bohr_per_angstrom;
        Atom atom(f, position);

        if (residueNameStr == "HOH" && parseOnlySolvent_) {
          structure.push_back(atom);
        }

        else if ((!parseOnlySolvent_ && residueNameStr != "HOH") && (overlayIdentifier == s || overlayIdentifier.empty())) {
          structure.push_back(atom);
        }
      }
      if (structure.size() == 0) {
        throw std::runtime_error("Error parsing the structure!");
      }
      structures.push_back(structure);
    }
  }
  return structures;
}

void PdbStreamHandler::write(std::ostream& os, const AtomCollection& atoms, const std::string& comment) {
  int N = atoms.size();
  os << comment << "\n";
  for (int i = 0; i < N; ++i) {
    std::string element = ElementInfo::symbol(atoms.getElement(i));
    auto position = atoms.getPosition(i) * Constants::angstrom_per_bohr;
    os << "ATOM" << std::setw(7) << std::right << i + 1 << "  " << std::setw(4) << std::left << element << std::setw(13)
       << std::left << "UNX" << std::setw(8) << std::right << std::fixed << std::setprecision(3) << position(0)
       << std::setw(8) << std::right << std::fixed << std::setprecision(3) << position(1) << std::setw(8) << std::right
       << std::fixed << std::setprecision(3) << position(2) << std::setw(7) << std::right << "    " << std::setw(6)
       << std::right << "     " << std::setw(11) << std::right << element << "\n";
  }
}

void PdbStreamHandler::setReadH(bool includeH) {
  includeH_ = includeH;
}

void PdbStreamHandler::parseOnlySolvent(bool parseOnlySolvent) {
  parseOnlySolvent_ = parseOnlySolvent;
}

void PdbStreamHandler::setSubstructureID(int substructureID) {
  substructureID_ = substructureID;
}

void PdbStreamHandler::extractOverlayIdentifiers(std::string line, PdbFileData& data) {
  std::string identifier = line.substr(16, 1);
  identifier.erase(std::remove(identifier.begin(), identifier.end(), ' '), identifier.end());
  if (!identifier.empty() && std::find(data.overlayIdentifiers.begin(), data.overlayIdentifiers.end(), identifier) ==
                                 data.overlayIdentifiers.end()) {
    data.overlayIdentifiers.push_back(identifier);
  }
}

bool PdbStreamHandler::isAtomLine(std::string line) {
  return (line.rfind("ATOM", 0) == 0) || (line.rfind("HETATM", 0) == 0);
}

bool PdbStreamHandler::isModelLine(std::string line) {
  return line.rfind("MODEL", 0) == 0;
}

void PdbStreamHandler::extractModels(std::istringstream& in, std::string& line, PdbFileData& data) {
  int modelNumber = 0;
  while (std::getline(in, line)) {
    if (isModelLine(line)) {
      data.nAtoms = 0;
      std::string atomBlock;
      while (!(line.rfind("ENDMDL", 0) == 0) && std::getline(in, line)) {
        // Extract atom block
        if (isAtomLine(line)) {
          data.nAtoms++;
          // Find overlaying substructures
          extractOverlayIdentifiers(line, data);
          std::istringstream atomLine(line);
          atomBlock.append(std::string(std::istreambuf_iterator<char>{atomLine}, {}));
          atomBlock.append("\n");
        }
      }
      data.atomBlocks.at(modelNumber) = atomBlock;
      modelNumber++;
    }
  }
}

void PdbStreamHandler::extractStructure(std::istringstream& in, std::string& line, PdbFileData& data) {
  while (std::getline(in, line)) {
    if (isAtomLine(line)) {
      // Find overlaying substructures
      extractOverlayIdentifiers(line, data);
      data.nAtoms++;
      std::istringstream atomLine(line);
      data.atomBlocks[0].append(std::string(std::istreambuf_iterator<char>{atomLine}, {}));
      data.atomBlocks[0].append("\n");
    }
  }
}

} // namespace Utils
} // namespace Scine
