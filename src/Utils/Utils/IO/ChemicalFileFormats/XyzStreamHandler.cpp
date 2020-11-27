/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/IO/ChemicalFileFormats/XyzStreamHandler.h"
#include "Utils/Bonds/BondOrderCollection.h"
#include "Utils/Constants.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Geometry/ElementInfo.h"
#include "Utils/Technical/ScopedLocale.h"
#include "boost/optional.hpp"
#include <iomanip>
#include <string>

namespace Scine {

namespace Utils {

// Helper to insert newlines
inline std::ostream& nl(std::ostream& os) {
  os << "\n";
  return os;
}

constexpr const char* XyzStreamHandler::model;

std::pair<Utils::AtomCollection, Utils::BondOrderCollection> XyzStreamHandler::read(std::istream& is,
                                                                                    const std::string& format) const {
  if (format != "xyz") {
    throw FormattedStreamHandler::FormatUnsupportedException();
  }

  return std::make_pair(read(is), Utils::BondOrderCollection());
}

void XyzStreamHandler::write(std::ostream& os, const std::string& format, const Utils::AtomCollection& atoms) const {
  if (format != "xyz") {
    throw FormattedStreamHandler::FormatUnsupportedException();
  }

  return write(os, atoms);
}

void XyzStreamHandler::write(std::ostream& /* os */, const std::string& format, const AtomCollection& /* atoms */,
                             const BondOrderCollection& /* bondOrders */
                             ) const {
  if (format != "xyz") {
    throw FormattedStreamHandler::FormatUnsupportedException();
  }

  throw FormattedStreamHandler::NoBondInformationException();
}

std::vector<XyzStreamHandler::FormatSupportPair> XyzStreamHandler::formats() const {
  return {{{"xyz", SupportType::ReadWrite}}};
}

bool XyzStreamHandler::formatSupported(const std::string& format, SupportType /* operation */
                                       ) const {
  if (format == "xyz") {
    return true;
  }

  return false;
}

std::string XyzStreamHandler::name() const {
  return XyzStreamHandler::model;
}

AtomCollection XyzStreamHandler::read(std::istream& is) {
  /* Make sure that the conversions are done in the C locale (otherwise, may
   * recognize commas as decimal separator).
   */
  is.imbue(std::locale("C"));

  ElementTypeCollection elements;
  PositionCollection positions;

  /* Extract the number of atoms. This number is principally enforced and
   * any deviation leads to a thrown exception
   * If the first line (without the ignored linebreak and spaces) is not a single integer
   * a FormatMismatch is thrown.
   * If the number of atoms in the positions after the read-in is different to the integer
   * a AtomNumberMismatch is thrown
   */
  int nAtoms;
  while (true) {
    std::string temp_str;
    std::getline(is, temp_str, is.widen('\n'));
    std::stringstream parser(temp_str);
    if (parser >> nAtoms && (parser >> std::ws).eof() && nAtoms >= 0)
      break; // success
    throw FormattedStreamHandler::FormatMismatchException();
  }

  // If this operation did not fail, we can preallocate the expected number of positions
  elements.reserve(nAtoms);
  positions.resize(nAtoms, Eigen::NoChange_t());

  // Skip the comment line
  is.ignore(std::numeric_limits<std::streamsize>::max(), is.widen('\n'));

  // Read atom lines
  std::string elementString;
  ElementType element;
  for (int i = 0; !is.eof(); ++i) {
    // Try to read an element string.
    is >> elementString;
    if (is.fail() && !is.eof()) {
      throw FormattedStreamHandler::FormatMismatchException();
    }
    else if (is.fail() && is.eof()) {
      break;
    }
    try {
      // Make first letter uppercase
      std::transform(std::begin(elementString), std::begin(elementString) + 1, std::begin(elementString), ::toupper);
      // Make other letters lowercase
      std::transform(std::begin(elementString) + 1, std::end(elementString), std::begin(elementString) + 1, ::tolower);
      element = ElementInfo::elementTypeForSymbol(elementString);
    }
    catch (...) {
      throw FormattedStreamHandler::FormatMismatchException();
    }
    elements.push_back(element);

    if (i < nAtoms) {
      // Read directly into the matrix
      is >> positions(i, 0) >> positions(i, 1) >> positions(i, 2);
    }
    else {
      // more atoms in file than specified at the top
      throw FormattedStreamHandler::AtomNumberMismatchException();
    }

    if (is.fail()) {
      throw FormattedStreamHandler::FormatMismatchException();
    }

    // Finish the line.
    is.ignore(std::numeric_limits<std::streamsize>::max(), is.widen('\n'));
  }

  if (elements.size() < nAtoms) {
    // fewer atoms in file than specified at the top
    throw FormattedStreamHandler::AtomNumberMismatchException();
  }

  /* Since all positions were read in as angstrom, but are internally
   * represented in bohr length units, we have to convert all of them:
   */
  positions *= Constants::bohr_per_angstrom;

  return {elements, positions};
}

void XyzStreamHandler::write(std::ostream& os, const AtomCollection& atoms) {
  os.imbue(std::locale("C"));

  // Write the atom count to the first line.
  os << std::setprecision(0) << std::fixed << atoms.size() << "\n";

  // Skip the comments line
  os << "\n";

  // Set floating point precision for positions
  os << std::setprecision(10);

  // Write positions
  const unsigned N = atoms.size();
  for (unsigned i = 0; i < N; ++i) {
    auto position = atoms.getPosition(i);

    os << std::left << std::setw(3) << ElementInfo::symbol(atoms.getElement(i));
    // clang-format off
    os << std::right
      << std::setw(16) << toAngstrom(Bohr(position.x()))
      << std::setw(16) << toAngstrom(Bohr(position.y()))
      << std::setw(16) << toAngstrom(Bohr(position.z()))
      << "\n";
    // clang-format on
  }
}

} // namespace Utils

} // namespace Scine
