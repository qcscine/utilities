/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "OrcaHessianOutputParser.h"
#include <Utils/ExternalQC/Exceptions.h>
#include <Utils/IO/Regex.h>
#include <fstream>
#include <iterator>
#include <regex>

namespace Scine {
namespace Utils {
namespace ExternalQC {

OrcaHessianOutputParser::OrcaHessianOutputParser(const std::string& filename) {
  extractContent(filename);
}

void OrcaHessianOutputParser::extractContent(const std::string& filename) {
  std::ifstream fin;
  fin.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  fin.open(filename);
  content_ = std::string(std::istreambuf_iterator<char>{fin}, {});
}

HessianMatrix OrcaHessianOutputParser::getHessian() const {
  std::istringstream iss(content_);

  readUntilHessianKeyword(iss);
  auto atomCount = getNumberAtomsFromHessianOutput(iss);
  Eigen::MatrixXd hessianMatrix(atomCount, atomCount);

  // The matrix is given five columns at a time, thus several blocks will be displayed after each other.
  // number blocks as function of number atoms: 0->0, 1->1, 5->1, 6->2
  auto numberBlocks = (atomCount + 4) / 5;

  for (int block = 0; block < numberBlocks; ++block) {
    ignoreFirstBlockLine(iss);
    int firstBlockColumnIndex = 5 * block;
    readOneBlock(iss, hessianMatrix, atomCount, firstBlockColumnIndex);
  }

  return hessianMatrix;
}

void OrcaHessianOutputParser::readUntilHessianKeyword(std::istream& in) const {
  std::string line;
  while (std::getline(in, line)) {
    if (line == "$hessian")
      return;
  }
  throw OutputFileParsingError("Could not find \"$hessian\" in hessian output file.");
}

int OrcaHessianOutputParser::getNumberAtomsFromHessianOutput(std::istream& in) const {
  std::string atomNumberLine;
  std::getline(in, atomNumberLine);
  std::regex regex(R"(^(\d+)$)");
  std::smatch matches;
  if (regex_search(atomNumberLine, matches, regex)) {
    return stoi(matches[1]);
  }
  throw OutputFileParsingError("Could not find the number of atoms in hessian output file.");
}

void OrcaHessianOutputParser::ignoreFirstBlockLine(std::istream& in) const {
  std::string firstBlockLine;
  std::getline(in, firstBlockLine);
}

void OrcaHessianOutputParser::readOneBlock(std::istream& in, Eigen::MatrixXd& m, int atomCount, int firstBlockColumnIndex) const {
  std::string line;
  for (int lineIndex = 0; lineIndex < atomCount; ++lineIndex) {
    std::getline(in, line);
    // Format:
    //   4    1.234567E-02   1.234567E-02   1.234567E-02   1.234567E-02   1.234567E-02
    std::string oneNumberRegex = R"(\s+)" + Regex::capturingFloatingPointNumber();

    std::string regexString = R"(\s*\d+)";
    int remainingAtoms = std::min(5, atomCount - firstBlockColumnIndex);

    for (int i = 0; i < remainingAtoms; ++i) {
      regexString += oneNumberRegex;
    }
    std::regex regex(regexString);
    std::smatch matches;
    if (regex_search(line, matches, regex)) {
      for (int i = 0; i < remainingAtoms; ++i) {
        m(lineIndex, firstBlockColumnIndex + i) = std::stod(matches[1 + i]);
      }
    }
    else {
      throw OutputFileParsingError("Error in parsing of Hessian, line: \"" + line + "\"");
    }
  }
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine