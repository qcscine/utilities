/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "TurbomoleHelper.h"
#include <Utils/ExternalQC/Exceptions.h>
#include <Utils/IO/NativeFilenames.h>
#include <boost/asio.hpp>
#include <boost/process.hpp>
#include <boost/process/async.hpp>
#include <regex>

namespace Scine {
namespace Utils {
namespace ExternalQC {

namespace bp = boost::process;
namespace bfs = boost::filesystem;

TurbomoleHelper::TurbomoleHelper(std::string& calculationDirectory, std::string& turbomoleExecutableBase)
  : calculationDirectory_(calculationDirectory), turbomoleExecutableBase_(turbomoleExecutableBase) {
}

bool TurbomoleHelper::jobWasSuccessful(std::istream& in) {
  std::regex r1("(ended normally)");
  std::smatch m1;

  // Get output
  std::string outputString;
  std::string line;
  while (std::getline(in, line))
    outputString += line;

  bool endedNormally = std::regex_search(outputString, m1, r1);
  return endedNormally;
}

void TurbomoleHelper::execute(std::string binaryName, bool outputToFile) {
  auto workingDirectory = bp::start_dir(calculationDirectory_);
  bp::ipstream err;
  std::string binary = NativeFilenames::combinePathSegments(turbomoleExecutableBase_, binaryName);

  if (outputToFile) {
    std::string outputFileName = binaryName + ".out";
    std::string outputFile = NativeFilenames::combinePathSegments(calculationDirectory_, outputFileName);
    // Delete output file before the calculation if it already exists.
    // Necessary since otherwise boost just pipes into the old output file which may lead to problems.
    bfs::remove(outputFile);
    bp::child c(binary, bp::std_out > outputFile, bp::std_err > err, workingDirectory);
    c.wait();
  }
  else {
    bp::ipstream out;
    bp::child c(binary, bp::std_out > out, bp::std_err > err, workingDirectory);
    c.wait(); // Wait for the process to exit
  }

  bool success = jobWasSuccessful(err);
  if (!success)
    throw std::runtime_error(binaryName + " calculation in Turbomole failed.");
}

void TurbomoleHelper::execute(std::string binaryName, std::string stdInFile) {
  auto workingDir = bp::start_dir(calculationDirectory_);
  bp::ipstream stderr;
  std::string binary = NativeFilenames::combinePathSegments(turbomoleExecutableBase_, binaryName);
  bp::child c(binary, bp::std_in<stdInFile, bp::std_out> bp::null, bp::std_err > stderr, workingDir);
  c.wait(); // Wait for the process to exit

  bool success = jobWasSuccessful(stderr);
  if (!success)
    throw std::runtime_error(binaryName + " session in Turbomole failed.");
}

void TurbomoleHelper::emptyFile(std::string file) {
  if (!file.empty()) {
    std::ofstream os;
    os.open(file, std::ofstream::out | std::ofstream::trunc);
    os.close();
  }
}

void TurbomoleHelper::mapBasisSetToTurbomoleStringRepresentation(std::string& basisSetString) {
  std::transform(std::begin(basisSetString), std::end(basisSetString), std::begin(basisSetString),
                 [](const auto c) { return std::tolower(c); });
  std::array<std::string, 4> basisSetIdentifiers = {"def2-", "def-", "cc-p", "aug-cc-p"};
  std::array<std::string, 3> basisSetToUpperCase = {"6-31g*", "sto-3g", "6-31g**"};

  // For a set of basis sets, partial transformation to uppercase is required
  bool isIdentified = false;
  std::string id;
  for (const auto& identifier : basisSetIdentifiers) {
    if (basisSetString.compare(0, identifier.length(), identifier) == 0) {
      isIdentified = true;
      id = identifier;
    }
  }

  if (isIdentified) {
    int length = id.length();
    std::string basisSetQuality = basisSetString.substr(basisSetString.find(id) + length);
    std::transform(std::begin(basisSetQuality), std::end(basisSetQuality), std::begin(basisSetQuality),
                   [](unsigned char u) { return std::toupper(u); });

    basisSetString = id + basisSetQuality;
  }
  // some other basis sets are uppercase only
  else if (std::find(std::begin(basisSetToUpperCase), std::end(basisSetToUpperCase), basisSetString) !=
           basisSetToUpperCase.end()) {
    std::transform(std::begin(basisSetString), std::end(basisSetString), std::begin(basisSetString),
                   [](unsigned char u) { return std::toupper(u); });
  }
  else
    throw std::runtime_error("Basis set " + basisSetString + " currently not supported by Turbomole calculator.");
}

void TurbomoleHelper::mapDftFunctionalToTurbomoleStringRepresentation(std::string& functionalString) {
  // make sure that functional is lowercase only
  std::transform(std::begin(functionalString), std::end(functionalString), std::begin(functionalString),
                 [](const auto c) { return std::tolower(c); });

  // check if functional should contain a delimiter
  std::unordered_map<std::string, std::string>::iterator it = correctedDftFunctionals_.find(functionalString);
  if (it != correctedDftFunctionals_.end()) {
    functionalString = it->second;
  }
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
