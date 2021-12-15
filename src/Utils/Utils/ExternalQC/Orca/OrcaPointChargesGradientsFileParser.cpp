/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "OrcaPointChargesGradientsFileParser.h"
#include <fstream>
#include <iterator>
#include <regex>

namespace Scine {
namespace Utils {
namespace ExternalQC {

OrcaPointChargesGradientsFileParser::OrcaPointChargesGradientsFileParser(const std::string& outputFileName) {
  extractContent(outputFileName);
}

void OrcaPointChargesGradientsFileParser::extractContent(const std::string& filename) {
  std::ifstream fin;
  fin.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  fin.open(filename);
  content_ = std::string(std::istreambuf_iterator<char>{fin}, {});
  fin.close();
}

GradientCollection OrcaPointChargesGradientsFileParser::getPointChargesGradients() const {
  std::istringstream iss(content_);
  std::string line;
  std::getline(iss, line);

  auto nRows = std::stoi(line);
  GradientCollection gradients(nRows, 3);
  std::getline(iss, line);

  for (int i = 0; i < nRows; ++i) {
    std::regex rgx("\\s+");
    std::sregex_token_iterator iter(line.begin(), line.end(), rgx, -1);
    std::sregex_token_iterator end;

    if (*iter == "") {
      iter++;
    }

    gradients(i, 0) = std::stod(*iter++);
    gradients(i, 1) = std::stod(*iter++);
    gradients(i, 2) = std::stod(*iter++);
    std::getline(iss, line);
  }

  return gradients;
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
