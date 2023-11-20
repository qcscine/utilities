/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/IO/MolecularTrajectoryIO.h"
#include "Utils/Constants.h"
#include "Utils/Geometry/ElementInfo.h"
#include "Utils/MolecularTrajectory.h"
#include <fstream>
#include <iomanip>

using namespace std;

namespace Scine {
namespace Utils {

void MolecularTrajectoryIO::write(format f, const std::string& fileName, const MolecularTrajectory& m) {
  ofstream fout;
  switch (f) {
    case format::binary:
      fout.open(fileName, ios_base::out | ios_base::trunc | ios_base::binary);
      break;
    case format::xyz:
      fout.open(fileName);
      break;
    case format::pdb:
      throw std::runtime_error("Only reading and NOT writing of pdb trajectories is supported at the moment");
  }
  if (!fout.is_open()) {
    throw std::runtime_error("Problem when opening/creating file " + fileName);
  }

  return write(f, fout, m);
}

void MolecularTrajectoryIO::write(format f, std::ostream& out, const MolecularTrajectory& m) {
  if (f == format::binary) {
    writeBinary(out, m);
  }
  else if (f == format::xyz) {
    writeXYZ(out, m);
  }
}

void MolecularTrajectoryIO::writeBinary(std::ostream& out, const MolecularTrajectory& m) {
  /*
   * Format : - number structures (int32_t)
              - number atoms (int32_t)
              - Z of each atom (int32_t)
              - coordinates of structures, in bohr (double)
   */

  int32_t nStructures = m.size();
  int32_t nAtoms = m.molecularSize();

  out.write(reinterpret_cast<char*>(&nStructures), sizeof(int32_t)); // NOLINT
  out.write(reinterpret_cast<char*>(&nAtoms), sizeof(int32_t));      // NOLINT

  // Atomic numbers
  const auto& elements = m.getElementTypes();
  for (int j = 0; j < nAtoms; ++j) {
    using ElementUnderlying = std::underlying_type<ElementType>::type;
    auto elementUnderlying = static_cast<ElementUnderlying>(elements[j]);
    out.write(reinterpret_cast<char*>(&elementUnderlying), sizeof(elementUnderlying)); // NOLINT
  }

  // Coordinates
  for (const auto& s : m) {
    for (int i = 0; i < s.rows(); i++) {
      out.write(reinterpret_cast<const char*>(s.row(i).data()), 3 * sizeof(s(i, 0))); // NOLINT
    }
  }
}

void MolecularTrajectoryIO::writeXYZ(std::ostream& out, const MolecularTrajectory& m) {
  out.imbue(std::locale("C"));
  const auto& elements = m.getElementTypes();
  bool hasEnergies = !m.getEnergies().empty();
  for (int i = 0; i < m.size(); ++i) {
    out << m.molecularSize() << endl;
    if (hasEnergies) {
      out << std::left << std::setw(10) << std::fixed << std::setprecision(10) << m.getEnergies()[i];
    }
    out << endl;
    for (int j = 0; j < m.molecularSize(); ++j) {
      writeXYZLine(out, elements[j], m[i].row(j));
    }
  }
}

void MolecularTrajectoryIO::writeXYZLine(std::ostream& out, ElementType e, const Position& p) {
  out << left << std::setw(3) << ElementInfo::symbol(e);
  out << right << setw(16) << fixed << setprecision(10) << p.x() * Constants::angstrom_per_bohr << right << setw(16)
      << fixed << setprecision(10) << p.y() * Constants::angstrom_per_bohr << right << setw(16) << fixed
      << setprecision(10) << p.z() * Constants::angstrom_per_bohr << endl;
}

MolecularTrajectory MolecularTrajectoryIO::read(MolecularTrajectoryIO::format f, const std::string& fileName) {
  ifstream fin;
  switch (f) {
    case format::binary:
      fin.open(fileName, ios_base::in | ios_base::binary);
      break;
    case format::xyz:
    case format::pdb:
      fin.open(fileName);
      break;
  }

  if (!fin.is_open()) {
    throw std::runtime_error("Problem when opening file " + fileName);
  }

  return read(f, fin);
}

MolecularTrajectory MolecularTrajectoryIO::read(MolecularTrajectoryIO::format f, std::istream& in) {
  // Use switch statement here so that the compiler complains if an option is not handled.
  switch (f) {
    case format::binary:
      return readBinary(in);
    case format::xyz:
      return readXYZ(in);
    case format::pdb:
      return readPdb(in);
  }
  throw std::runtime_error("Unsupported format to read MolecularTrajectory from");
}

MolecularTrajectory MolecularTrajectoryIO::readBinary(std::istream& in) {
  /*
   * Format : - number structures (int32_t)
              - number atoms (int32_t)
              - Z of each atom (int32_t)
              - coordinates of structures, in bohr (double)
   */
  int32_t nStructures;
  int32_t nAtoms;
  in.read(reinterpret_cast<char*>(&nStructures), sizeof(int32_t)); // NOLINT
  in.read(reinterpret_cast<char*>(&nAtoms), sizeof(int32_t));      // NOLINT

  ElementTypeCollection elements(nAtoms);

  // Atomic numbers
  for (int j = 0; j < nAtoms; ++j) {
    using ElementUnderlying = std::underlying_type<ElementType>::type;
    ElementUnderlying elementUnderlying;
    in.read(reinterpret_cast<char*>(&elementUnderlying), sizeof(elementUnderlying)); // NOLINT
    elements[j] = static_cast<ElementType>(elementUnderlying);
  }

  MolecularTrajectory m;
  m.setElementTypes(elements);
  m.resize(nStructures);

  // Coordinates
  PositionCollection tempPositions(nAtoms, 3);
  for (auto& s : m) {
    for (int j = 0; j < nAtoms; ++j) {
      in.read(reinterpret_cast<char*>(s.row(j).data()), 3 * sizeof(s(j, 0))); // NOLINT
    }
  }

  return m;
}

MolecularTrajectory MolecularTrajectoryIO::readXYZ(std::istream& in) {
  in.imbue(std::locale("C"));
  bool firstDone = false;
  ElementTypeCollection elements;
  MolecularTrajectory m;
  std::vector<double> energies;
  while (!in.eof()) {
    int nAtoms;
    in >> nAtoms;
    if (in.eof()) {
      break;
    }

    std::string unnecessaryLine;
    getline(in, unnecessaryLine);
    std::string commentLine;
    getline(in, commentLine);
    bool gotEnergy = !commentLine.empty() && commentLine.find_first_not_of("-.0123456789") == string::npos;
    if (gotEnergy) {
      // in case we got funky comment, that really isn't a double
      try {
        energies.push_back(std::stod(commentLine));
      }
      catch (...) {
        ;
      }
    }
    std::string element;
    PositionCollection s(nAtoms, 3);
    for (int i = 0; i < nAtoms; ++i) {
      in >> element;

      // Add elements only the first time
      if (!firstDone) {
        // Strip suffix numbers
        element.erase(
            std::remove_if(element.begin(), element.end(), [](std::string::value_type ch) { return isalpha(ch) == 0; }),
            element.end());
        elements.push_back(ElementInfo::elementTypeForSymbol(element));
      }
      double x, y, z;
      in >> x >> y >> z;
      s(i, 0) = x;
      s(i, 1) = y;
      s(i, 2) = z;
      getline(in, unnecessaryLine);
    }
    firstDone = true;
    m.push_back(s * Constants::bohr_per_angstrom);
  }
  m.setElementTypes(elements);
  if (energies.size() == static_cast<unsigned>(m.size())) {
    m.setEnergies(energies);
  }

  return m;
}

MolecularTrajectory MolecularTrajectoryIO::readPdb(std::istream& in) {
  MolecularTrajectory trajectory;
  ElementTypeCollection elements;

  /*
   * Example pdb file:
   * REMARK   1 CREATED WITH OPENMM 7.7, 2023-08-11
   * CRYST1   48.751   43.841   46.944  90.00  90.00  90.00 P 1           1
   * MODEL        1
   * HETATM    1  C    R0 A   1      21.588  21.544  23.151  1.00  0.00           C
   * HETATM    2  C    R1 A   2      23.127  21.575  23.151  1.00  0.00           C
   * ...
   * HETATM 2518  O    R6 A2518       0.069  40.026  44.070  1.00  0.00           O
   * HETATM 2519  H    R7 A2519      17.352  19.012  11.461  1.00  0.00           H
   * HETATM 2520  H    R7 A2520      40.264  31.406   4.836  1.00  0.00           H
   * TER    2521       R7 A2520
   * ENDMDL
   * MODEL        2
   * HETATM    1  C    R0 A   1      21.651  21.607  23.218  1.00  0.00           C
   * HETATM    2  C    R1 A   2      23.194  21.638  23.218  1.00  0.00           C
   * HETATM    3  O    R2 A   3      23.935  20.444  23.218  1.00  0.00           O
   * ...
   * TER    2521       R7 A2520
   * ENDMMDL
   * END
   */

  std::string skip;
  std::string coordinateLine;
  // Skip first two lines.
  getline(in, skip);
  getline(in, skip);
  unsigned int nAtoms = 0;
  bool firstPass = true;
  while (!in.eof()) {
    getline(in, coordinateLine);
    if (coordinateLine.find("END") != std::string::npos) {
      break;
    }
    std::vector<Eigen::Vector3d> coords;
    while (!in.eof()) {
      getline(in, coordinateLine);
      if (coordinateLine.find("END") != std::string::npos) {
        break;
      }
      if (coordinateLine.find("MODEL") != std::string::npos || coordinateLine.find("TER") != std::string::npos ||
          coordinateLine.find("CONNECT") != std::string::npos) {
        continue;
      }
      const std::string atomLabel = coordinateLine.substr(0, 6);
      if (atomLabel != "ATOM" && atomLabel != "HETATM") {
        throw runtime_error("Unexpected label " + atomLabel + " for atom. Only ATOM and HETATM are allowed in pdb files!");
      }
      if (firstPass) {
        std::string element = coordinateLine.substr(76, 2);
        // Remove white space if any are there.
        const auto spacePosition = element.find(' ');
        if (spacePosition != std::string::npos) {
          element = coordinateLine.substr(77, 1);
        }
        elements.push_back(ElementInfo::elementTypeForSymbol(element));
      }
      const double xCoordinate = std::stod(coordinateLine.substr(30, 8));
      const double yCoordinate = std::stod(coordinateLine.substr(38, 8));
      const double zCoordinate = std::stod(coordinateLine.substr(46, 8));
      coords.emplace_back(xCoordinate, yCoordinate, zCoordinate);
    }
    if (firstPass) {
      trajectory.setElementTypes(elements);
      firstPass = false;
      nAtoms = coords.size();
    }
    if (nAtoms != coords.size()) {
      throw runtime_error("Inconsistent number of atoms in trajectory file.");
    }
    PositionCollection positions(coords.size(), 3);
    for (unsigned int iCoord = 0; iCoord < coords.size(); ++iCoord) {
      positions.row(iCoord) = coords[iCoord];
    }
    trajectory.push_back(positions);
  }
  return trajectory;
}

} /* namespace Utils */
} /* namespace Scine */
