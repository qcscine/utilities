/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/IO/MolecularTrajectoryIO.h"
#include "Utils/Constants.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Geometry/ElementInfo.h"
#include "Utils/MolecularTrajectory.h"
#include <fstream>
#include <iomanip>

using namespace std;

namespace Scine {
namespace Utils {

void MolecularTrajectoryIO::write(format f, const std::string& fileName, const MolecularTrajectory& m) {
  ofstream fout;

  if (f == format::binary)
    fout.open(fileName, ios_base::out | ios_base::trunc | ios_base::binary);
  else
    fout.open(fileName);

  if (!fout.is_open())
    throw std::runtime_error("Problem when opening/creating file " + fileName);

  return write(f, fout, m);
}

void MolecularTrajectoryIO::write(format f, std::ostream& out, const MolecularTrajectory& m) {
  if (f == format::binary)
    writeBinary(out, m);
  else if (f == format::xyz)
    writeXYZ(out, m);
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
  for (auto& s : m) {
    for (int i = 0; i < s.rows(); i++) {
      out.write(reinterpret_cast<const char*>(s.row(i).data()), 3 * sizeof(s(i, 0))); // NOLINT
    }
  }
}

void MolecularTrajectoryIO::writeXYZ(std::ostream& out, const MolecularTrajectory& m) {
  out.imbue(std::locale("C"));
  const auto& elements = m.getElementTypes();
  for (int i = 0; i < m.size(); ++i) {
    out << m.molecularSize() << endl << endl;
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

  if (f == format::binary)
    fin.open(fileName, ios_base::in | ios_base::binary);
  else
    fin.open(fileName);

  if (!fin.is_open()) {
    throw std::runtime_error("Problem when opening file " + fileName);
  }

  return read(f, fin);
}

MolecularTrajectory MolecularTrajectoryIO::read(MolecularTrajectoryIO::format f, std::istream& in) {
  if (f == format::binary)
    return readBinary(in);
  if (f == format::xyz)
    return readXYZ(in);

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
  while (!in.eof()) {
    int nAtoms;
    in >> nAtoms;
    if (in.eof())
      break;

    std::string unnecessaryLine;
    getline(in, unnecessaryLine);
    getline(in, unnecessaryLine);
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

  return m;
}

} /* namespace Utils */
} /* namespace Scine */
