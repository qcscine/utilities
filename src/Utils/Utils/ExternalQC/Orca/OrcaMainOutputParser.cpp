/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "OrcaMainOutputParser.h"
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

int OrcaMainOutputParser::getNumberAtoms() const {
  std::regex regex(R"(Number of atoms +\.+ +(\d+))");
  std::smatch matches;
  bool b = std::regex_search(content_, matches, regex);
  if (!b)
    throw OutputFileParsingError("Number of atoms could not be read from ORCA output");
  return std::stoi(matches[1]);
}

GradientCollection OrcaMainOutputParser::getGradients() const {
  int nAtoms = getNumberAtoms();

  // first go to section about gradients
  std::regex r1(R"((CARTESIAN GRADIENT|The cartesian gradient:))");
  std::smatch m1;
  bool b = std::regex_search(content_, m1, r1);
  if (!b)
    throw OutputFileParsingError("Gradient section could not be found in ORCA output");
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
      throw OutputFileParsingError("Gradient could not be found in ORCA output");
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
    throw OutputFileParsingError("Energy could not be read from ORCA output");
  return std::stod(matches[1]);
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
