/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "TurbomoleMainOutputParser.h"
#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/ExternalQC/Exceptions.h>
#include <Utils/IO/Regex.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <boost/process.hpp>
#include <fstream>
#include <iostream>
#include <iterator>
#include <regex>

namespace Scine {
namespace Utils {
namespace ExternalQC {

TurbomoleMainOutputParser::TurbomoleMainOutputParser(TurbomoleFiles& files) : files_(files) {
  extractContent(files.outputFile);
}

void TurbomoleMainOutputParser::extractContent(const std::string& filename) {
  std::ifstream fin;
  fin.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  fin.open(filename);
  content_ = std::string(std::istreambuf_iterator<char>{fin}, {});
  fin.close();
}

int TurbomoleMainOutputParser::getNumberAtoms() const {
  std::ifstream atoms;
  atoms.open(files_.coordFile);
  std::string line;
  int nAtoms = 0;

  std::regex r1(R"(([-\.0-9])+ +([-\.0-9])+ +([-\.0-9])+ +([A-Za-z])+)");
  std::smatch m1;

  while (std::getline(atoms, line)) {
    if (std::regex_search(line, m1, r1))
      nAtoms++;
  }
  atoms.close();
  return nAtoms;
}

int TurbomoleMainOutputParser::getNumberOfNonZeroPointCharges() const {
  std::ifstream pc;
  pc.open(files_.pointChargesFile);
  std::string line;

  int numPointCharges = 0;
  while (std::getline(pc, line)) {
    std::vector<std::string> lineSplitted;
    boost::split(lineSplitted, line, boost::is_any_of(" "), boost::token_compress_on);
    try {
      std::stod(lineSplitted[0]);
      std::stod(lineSplitted[1]);
      std::stod(lineSplitted[2]);
      double charge = std::stod(lineSplitted[3]);
      if (charge != 0.0) {
        numPointCharges++;
      }
    }
    catch (std::exception& e) {
      throw std::runtime_error("Point charges file " + files_.pointChargesFile + " has incorrect format!");
    }
  }
  pc.close();
  return numPointCharges;
}

GradientCollection TurbomoleMainOutputParser::getGradients() const {
  int nAtoms = getNumberAtoms();
  GradientCollection gc;
  gc.resize(nAtoms, 3);

  std::ifstream in(files_.gradientFile);

  for (int i = 0; i < nAtoms + 2; i++) {
    std::string line;
    std::getline(in, line);
  }

  for (int i = 0; i < nAtoms; i++) {
    std::array<std::string, 3> s;
    in >> s[0] >> s[1] >> s[2];
    for (std::string& j : s) {
      const unsigned long e = j.find_first_of("Dd");
      if (e != std::string::npos) {
        j[e] = 'E';
      }
    }
    try {
      gc(i, 0) = std::stod(s[0]);
      gc(i, 1) = std::stod(s[1]);
      gc(i, 2) = std::stod(s[2]);
    }
    catch (...) {
      throw OutputFileParsingError("Gradient file " + files_.gradientFile + " could not be parsed.");
    }
  }

  return gc;
}

GradientCollection TurbomoleMainOutputParser::getPointChargesGradients() const {
  int nPointCharges = getNumberOfNonZeroPointCharges();
  if (nPointCharges == 0) {
    throw std::runtime_error("Error parsing the point charges!");
  }
  GradientCollection gc;
  gc.resize(nPointCharges, 3);

  std::ifstream in(files_.pointChargeGradientFile);
  // skip first line
  std::string line;
  std::getline(in, line);
  for (int i = 0; i < nPointCharges; i++) {
    std::array<std::string, 3> s;
    in >> s[0] >> s[1] >> s[2];
    for (std::string& j : s) {
      const unsigned long e = j.find_first_of("Dd");
      if (e != std::string::npos) {
        j[e] = 'E';
      }
    }
    try {
      gc(i, 0) = std::stod(s[0]);
      gc(i, 1) = std::stod(s[1]);
      gc(i, 2) = std::stod(s[2]);
    }
    catch (...) {
      throw OutputFileParsingError("Gradient file " + files_.pointChargeGradientFile + " for n point charges " +
                                   std::to_string(nPointCharges) + " could not be parsed due to line " + s[0] + s[1] +
                                   s[2] + " at line " + std::to_string(i));
    }
  }

  return gc;
}

double TurbomoleMainOutputParser::getEnergy() const {
  std::ifstream in;
  in.open(files_.energyFile);
  std::string content = std::string(std::istreambuf_iterator<char>{in}, {});
  in.close();

  std::regex energyRegex(R"(\d+ +([-\.0-9]+) +[-\.0-9]+ +[-\.0-9]+)");

  std::string line;
  std::smatch m;
  bool found = std::regex_search(content, m, energyRegex);
  if (!found)
    throw OutputFileParsingError("Energy could not be read from Turbomole output.");

  return std::stod(m[1]);
}

double TurbomoleMainOutputParser::getExcitedStateEnergy(int state) const {
  std::ifstream in;
  in.open(files_.escfFile);
  std::string content = std::string(std::istreambuf_iterator<char>{in}, {});
  in.close();

  std::regex regex(std::string("\\s+") + std::to_string(state) + " a excitation\\s+Total energy:\\s+(-?)\\d+\\.\\d+");
  std::smatch m;
  bool found = std::regex_search(content, m, regex);
  if (!found) {
    throw OutputFileParsingError("Excited state energy could not be read from Turbomole output.");
  }

  std::string line = m[0].str();
  std::string delimiter = ":";
  int index = line.find(delimiter) + delimiter.size();
  std::string substring = line.substr(index, line.size());

  return std::stod(substring);
}

Utils::BondOrderCollection TurbomoleMainOutputParser::getBondOrders() const {
  int nAtoms = getNumberAtoms();
  BondOrderCollection bondOrders(nAtoms);

  // First go to section about Hirshfeld charges
  std::regex r1(R"(WIBERG BOND INDICES)");
  std::smatch m1;
  bool b = std::regex_search(content_, m1, r1);
  if (!b) {
    throw OutputFileParsingError("Bond Orders section could not be found in TURBOMOLE output.");
  }

  /* Format:
   * [ElementType1] [index1] [ElementType2] [index2] [bondOrder]
   */

  std::regex bondOrderRegex(R"(([A-Za-z]+) +([0-9]+) +-+ +([A-Za-z]+) +([0-9]+) +([-\.0-9]+))");
  std::smatch m2;
  std::istringstream contentStream(content_);
  std::string line;
  while (std::getline(contentStream, line)) {
    if (std::regex_search(line, m2, bondOrderRegex)) {
      int i = std::stoi(m2[2]);
      int j = std::stoi(m2[4]);
      double bo = std::stod(m2[5]);
      bondOrders.setOrder(i - 1, j - 1, bo);
    }
  }

  return bondOrders;
}

void TurbomoleMainOutputParser::checkForErrors() const {
  std::regex regex("(convergence criteria cannot be satisfied)");
  std::smatch match;
  if (std::regex_search(content_, match, regex)) {
    throw ScfNotConvergedError("SCF in Turbomole calculation did not converge.");
  }
}

std::vector<double> TurbomoleMainOutputParser::getLoewdinCharges() const {
  std::vector<double> charges;
  int nAtoms = getNumberAtoms();

  // First go to section about Hirshfeld charges
  std::regex r1(R"(atom      charge)");
  std::smatch m1;
  bool b = std::regex_search(content_, m1, r1);
  if (!b) {
    throw OutputFileParsingError("Loewdin charges section could not be found in TURBOMOLE output.");
  }
  // Set iterator where to start search for Hirshfeld charges
  auto it = m1[0].second;

  /* Format:
   * [ElementType] [LoewdinCharge] [uselessFloatingPointNumbers] ..
   */
  std::regex r2(R"([A-Za-z]+ +([-\.0-9]+) +[-\.0-9]+)");
  std::smatch m2;
  for (int i = 0; i < nAtoms; ++i) {
    if (std::regex_search(it, content_.end(), m2, r2)) {
      auto charge = std::stod(m2[1]);
      charges.push_back(charge);
      it = m2[0].second;
    }
    else {
      throw OutputFileParsingError("Loewdin charges could not be found in TURBOMOLE output.");
    }
  }
  return charges;
}

HessianMatrix TurbomoleMainOutputParser::getHessian() const {
  int nAtoms = getNumberAtoms();
  Eigen::MatrixXd hessianMatrix(3 * nAtoms, 3 * nAtoms);

  std::ifstream in;
  in.open(files_.hessianFile);
  std::string line;
  std::vector<double> vals;
  while (std::getline(in, line)) {
    if ((line.find("$hessian") != std::string::npos) || (line.find("$end") != std::string::npos))
      continue;
    else {
      std::istringstream iss(line);
      std::vector<std::string> result{std::istream_iterator<std::string>(iss), {}};
      for (auto s : result) {
        if (s.end() == std::find_if(s.begin(), s.end(), [](unsigned char c) -> bool { return !isdigit(c); }))
          continue;
        else
          vals.push_back(std::stod(s));
      }
    }
  }

  in.close();

  Eigen::MatrixXd matrix = Eigen::Map<Eigen::MatrixXd>(vals.data(), 3 * nAtoms, 3 * nAtoms);
  if (!(matrix.rows() == 3 * nAtoms) || !(matrix.cols() == 3 * nAtoms))
    throw std::runtime_error("Parsing of the hessian matrix failed.");
  if (!matrix.isApprox(matrix.transpose()))
    throw std::runtime_error("Parsed hessian is not symmetric.");

  return matrix;
}

int TurbomoleMainOutputParser::getSymmetryNumber() const {
  // TODO: Turbomole does not print the symmetry number
  return 1;
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
