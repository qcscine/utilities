/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "OrcaMainOutputParser.h"
#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/ExternalQC/Exceptions.h>
#include <Utils/IO/Regex.h>
#include <fstream>
#include <iterator>
#include <regex>

namespace Scine {
namespace Utils {
namespace ExternalQC {

OrcaMainOutputParser::OrcaMainOutputParser(const std::string& filename) {
  extractContent(filename);
}

void OrcaMainOutputParser::extractContent(const std::string& filename) {
  std::ifstream fin;
  fin.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  fin.open(filename);
  content_ = std::string(std::istreambuf_iterator<char>{fin}, {});
}

void OrcaMainOutputParser::checkForErrors() const {
  std::regex r1(R"(E\s?R\s?R\s?O\s?R\s?)");
  std::smatch m1;
  if (std::regex_search(content_, m1, r1))
    throw OutputFileParsingError("ORCA encountered an error during the calculation.");
}

int OrcaMainOutputParser::getNumberAtoms() const {
  int nAtoms = 0;

  std::istringstream contentStream(content_);
  std::string line;
  bool coordinatesStarted = false;
  bool coordinatesEnded = false;
  while (std::getline(contentStream, line)) {
    if (!coordinatesStarted) {
      if (line.find("CARTESIAN COORDINATES (ANGSTROEM)") != std::string::npos) {
        coordinatesStarted = true;
        // Do not try to parse this line
        continue;
      }
    }

    if (coordinatesStarted && !coordinatesEnded) {
      if (line.empty()) {
        // This indicates the end of the coordinates block
        coordinatesEnded = true;
        continue;
      }
      nAtoms++;
    }
  }

  if (!coordinatesStarted)
    throw OutputFileParsingError("Number of atoms could not be read from ORCA output.");

  return nAtoms - 1; // Removing one for a wrongly counted separator line
}

GradientCollection OrcaMainOutputParser::getGradients() const {
  int nAtoms = getNumberAtoms();

  // first go to section about gradients
  std::regex r1(R"((CARTESIAN GRADIENT|The cartesian gradient:))");
  std::smatch m1;
  bool b = std::regex_search(content_, m1, r1);
  if (!b)
    throw OutputFileParsingError("Gradient section could not be found in ORCA output.");
  // set iterator where to start search for gradients
  auto it = m1[0].second;

  GradientCollection gc;
  gc.resize(nAtoms, 3);
  std::string spacesAndNumber = " +" + Regex::capturingFloatingPointNumber();
  std::string regexString = R"(\d+ +\S+ +:)" + spacesAndNumber + spacesAndNumber + spacesAndNumber;
  std::regex r2(regexString);
  std::smatch m2;
  for (int i = 0; i < nAtoms; ++i) {
    if (std::regex_search(it, content_.end(), m2, r2)) {
      double x = std::stod(m2[1]);
      double y = std::stod(m2[2]);
      double z = std::stod(m2[3]);
      gc.row(i) = Gradient(x, y, z);
      it = m2[0].second;
    }
    else {
      throw OutputFileParsingError("Gradient could not be found in ORCA output.");
    }
  }
  return gc;
}

double OrcaMainOutputParser::getEnergy() const {
  std::string regexString = "FINAL SINGLE POINT ENERGY +" + Regex::capturingFloatingPointNumber();
  std::regex regex(regexString);
  std::smatch matches;
  bool b = std::regex_search(content_, matches, regex);
  if (!b)
    throw OutputFileParsingError("Energy could not be read from ORCA output.");
  return std::stod(matches[1]);
}

Utils::BondOrderCollection OrcaMainOutputParser::getBondOrders() const {
  int nAtoms = getNumberAtoms();
  Utils::BondOrderCollection bondOrders(nAtoms);

  const std::regex bondOrderRegex(R"(\w\(\s*(-?\d+)-\w+\s*,\s*(-?\d+)-\w+\s*\)\s+:\s+(-?\d+\.\d+)\s+)");

  std::istringstream contentStream(content_);
  std::string line;
  bool parseLine = false;
  bool parsedBondOrders = false;
  while (std::getline(contentStream, line)) {
    if (line.find("Mayer bond orders larger than") != std::string::npos) {
      // Ensure that last Mayer bond order block is stored
      if (parsedBondOrders) {
        bondOrders.setZero();
      }
      parseLine = true;
      // Do not try to parse this line, the next is the first interesting one
      continue;
    }

    if (parseLine) {
      if (line.empty()) {
        // This indicates the end of the Mayer bond order block
        parseLine = false;
        parsedBondOrders = true;
        continue;
      }

      // Parse the line with a regex
      auto line_begin = std::sregex_iterator(std::begin(line), std::end(line), bondOrderRegex);
      auto line_end = std::sregex_iterator();

      for (auto it = line_begin; it != line_end; ++it) {
        // Extract the capture groups. Stoi and stod all throw on error
        int i = std::stoi(it->str(1));
        int j = std::stoi(it->str(2));
        double bo = std::stod(it->str(3));

        bondOrders.setOrder(i, j, bo);
      }
    }
  }

  if (!parsedBondOrders)
    throw OutputFileParsingError("Bond orders could not be read from ORCA output.");
  return bondOrders;
}

std::vector<double> OrcaMainOutputParser::getHirshfeldCharges() const {
  int nAtoms = getNumberAtoms();
  std::vector<double> charges;

  // First go to section about Hirshfeld charges
  std::regex r1(R"(HIRSHFELD +ANALYSIS)");
  std::smatch m1;
  bool b = std::regex_search(content_, m1, r1);
  if (!b)
    throw OutputFileParsingError("Hirshfeld charges section could not be found in ORCA output.");
  // Set iterator where to start search for Hirshfeld charges
  auto it = m1[0].second;

  /* Format:
   * [atom index] [ElementType] [HirshfeldCharge] [uselessFloatingPointNumber]
   */
  std::regex r2(R"(\d+ +[A-Za-z]+ +([-\.0-9]+) +[-\.0-9]+)");
  std::smatch m2;
  for (int i = 0; i < nAtoms; ++i) {
    if (std::regex_search(it, content_.end(), m2, r2)) {
      double charge = std::stod(m2[1]);
      charges.push_back(charge);
      it = m2[0].second;
    }
    else {
      throw OutputFileParsingError("Hirshfeld charges could not be found in ORCA output.");
    }
  }

  return charges;
}

double OrcaMainOutputParser::getTemperature() const {
  std::string regexString = "Temperature+\\s+...\\s+" + Regex::capturingFloatingPointNumber();
  std::regex regex(regexString);
  std::smatch matches;
  bool b = std::regex_search(content_, matches, regex);
  if (!b)
    throw OutputFileParsingError("Temperature could not be read from ORCA output.");
  return std::stod(matches[1]);
}

double OrcaMainOutputParser::getEnthalpy() const {
  std::string regexString = "Total enthalpy+\\s+...\\s+" + Regex::capturingFloatingPointNumber();
  std::regex regex(regexString);
  std::smatch matches;
  bool b = std::regex_search(content_, matches, regex);
  if (!b)
    throw OutputFileParsingError("Enthalpy could not be read from ORCA output.");
  return std::stod(matches[1]);
}

double OrcaMainOutputParser::getEntropy() const {
  std::string regexString = "Total entropy correction+\\s+...\\s+" + Regex::capturingFloatingPointNumber();
  std::regex regex(regexString);
  std::smatch matches;
  bool b = std::regex_search(content_, matches, regex);
  if (!b)
    throw OutputFileParsingError("Entropy could not be read from ORCA output.");
  return -std::stod(matches[1]) / getTemperature();
}

double OrcaMainOutputParser::getZeroPointVibrationalEnergy() const {
  std::string regexString = "Non-thermal \\(ZPE\\) correction+\\s+...\\s+" + Regex::capturingFloatingPointNumber();
  std::regex regex(regexString);
  std::smatch matches;
  bool b = std::regex_search(content_, matches, regex);
  if (!b)
    throw OutputFileParsingError("Zero point vibrational energy could not be read from ORCA output.");
  return std::stod(matches[1]);
}

double OrcaMainOutputParser::getGibbsFreeEnergy() const {
  // In Orca 4.1.0 the employed word is 'enthalpy' and in Orca 4.2.0 it is 'energy'
  std::string regexString = "Final Gibbs free (?:enthalpy|energy)+\\s+...\\s+" + Regex::capturingFloatingPointNumber();
  std::regex regex(regexString);
  std::smatch matches;
  bool b = std::regex_search(content_, matches, regex);
  if (!b)
    throw OutputFileParsingError("Gibbs free energy could not be read from ORCA output.");
  return std::stod(matches[1]);
}

double OrcaMainOutputParser::getSymmetryNumber() const {
  std::string regexString = "Point Group:\\s+[a-zA-Z0-9]*\\s*,\\s+Symmetry Number:\\s+" + Regex::capturingIntegerNumber();
  std::regex regex(regexString);
  std::smatch matches;
  bool b = std::regex_search(content_, matches, regex);
  if (!b)
    throw OutputFileParsingError("Symmetry number could not be read from ORCA output.");
  return std::stod(matches[1]);
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
