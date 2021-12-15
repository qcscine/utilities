/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "GaussianOutputParser.h"
#include <Utils/ExternalQC/Exceptions.h>
#include <Utils/IO/Regex.h>
#include <fstream>
#include <iterator>
#include <regex>

namespace Scine {
namespace Utils {
namespace ExternalQC {

GaussianOutputParser::GaussianOutputParser(const std::string& filename) {
  extractContent(filename);
}

void GaussianOutputParser::extractContent(const std::string& filename) {
  std::ifstream fin;
  fin.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  fin.open(filename);
  content_ = std::string(std::istreambuf_iterator<char>{fin}, {});
  fin.close();
}

int GaussianOutputParser::getNumberAtoms() const {
  std::regex regex(R"(NAtoms= +(\d+))");
  std::smatch matches;
  bool b = std::regex_search(content_, matches, regex);
  if (!b) {
    throw OutputFileParsingError("Number of atoms could not be read from Gaussian output");
  }
  return std::stoi(matches[1]);
}

GradientCollection GaussianOutputParser::getGradients() const {
  int nAtoms = getNumberAtoms();

  // first go to section about forces
  std::regex r1(R"(Center +Atomic +Forces \(Hartrees/Bohr\))");
  std::smatch m1;
  bool b = std::regex_search(content_, m1, r1);
  if (!b) {
    throw OutputFileParsingError("Force section could not be found in GAUSSIAN output");
  }
  // set iterator where to start search for gradients
  auto it = m1[0].second;

  Utils::GradientCollection gc;
  gc.resize(nAtoms, 3);
  /* Format:
   * [atom index] [Z] [Fx] [Fy] [Fz]
   */
  std::regex r2(R"(\d+ +\d+ +([-\.0-9]+) +([-\.0-9]+) +([-\.0-9]+))");
  std::smatch m2;
  for (int i = 0; i < nAtoms; ++i) {
    if (std::regex_search(it, content_.end(), m2, r2)) {
      double x = std::stod(m2[1]);
      double y = std::stod(m2[2]);
      double z = std::stod(m2[3]);
      gc.row(i) = Eigen::RowVector3d(-x, -y, -z);
      it = m2[0].second;
    }
    else {
      throw OutputFileParsingError("Forces could not be found in GAUSSIAN output");
    }
  }
  return gc;
}

double GaussianOutputParser::getEnergy() const {
  std::regex regex(R"(SCF Done:  E\(.+\) = +([-\.0-9]+) +[Aa]\.[Uu]\.)");
  std::smatch matches;
  bool b = std::regex_search(content_, matches, regex);
  if (!b) {
    throw OutputFileParsingError("Energy could not be read from GAUSSIAN output");
  }
  return std::stod(matches[1]);
}

std::vector<double> GaussianOutputParser::getCM5Charges() const {
  int nAtoms = getNumberAtoms();
  std::vector<double> cm5Charges;

  // first go to section about cm5 charges
  std::regex r1(R"(Hirshfeld +charges, +spin +densities, +dipoles, +and +CM5 +charges)");
  std::smatch m1;
  bool b = std::regex_search(content_, m1, r1);
  if (!b) {
    throw OutputFileParsingError("CM5 charges section could not be found in GAUSSIAN output.");
  }
  // set iterator where to start search for cm5 charges
  auto it = m1[0].second;

  /* Format:
   * [atom index] [ElementType] [uselessCharge] [uselessCharge] [uselessCharge] [uselessCharge] [uselessCharge] [cm5
   * charge]
   */
  std::regex r2(R"(\d+ +[A-Za-z]+ +[-\.0-9]+ +[-\.0-9]+ +[-\.0-9]+ +[-\.0-9]+ +[-\.0-9]+ +([-\.0-9]+))");
  std::smatch m2;
  for (int i = 0; i < nAtoms; ++i) {
    if (std::regex_search(it, content_.end(), m2, r2)) {
      double charge = std::stod(m2[1]);
      cm5Charges.push_back(charge);
      it = m2[0].second;
    }
    else {
      throw OutputFileParsingError("CM5 charges could not be found in GAUSSIAN output.");
    }
  }

  return cm5Charges;
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
