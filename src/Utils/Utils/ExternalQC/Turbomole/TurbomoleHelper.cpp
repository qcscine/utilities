/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "TurbomoleHelper.h"
#include "Utils/MSVCCompatibility.h"
#include <Utils/ExternalQC/Exceptions.h>
#include <Utils/IO/NativeFilenames.h>
#include <boost/asio.hpp>
#include <boost/process.hpp>
#include <boost/process/async.hpp>
#include <fstream>
#include <regex>

namespace Scine {
namespace Utils {
namespace ExternalQC {

namespace bp = boost::process;
namespace bfs = boost::filesystem;

TurbomoleHelper::TurbomoleHelper(std::string& calculationDirectory, std::string& turbomoleExecutableBase)
  : calculationDirectory_(calculationDirectory), turbomoleExecutableBase_(turbomoleExecutableBase) {
}

bool TurbomoleHelper::jobWasSuccessful(std::istream& out, std::string phrase) {
  std::regex r1(phrase);
  std::smatch m1;

  // Get output
  std::string outputString;
  std::string line;
  while (std::getline(out, line))
    outputString += line;
  bool endedNormally = std::regex_search(outputString, m1, r1);
  return endedNormally;
}

void TurbomoleHelper::execute(std::string binaryName, bool outputToFile, std::string outputFileName) {
  auto workingDirectory = bp::start_dir(calculationDirectory_);
  bp::ipstream err;
  std::string binary = NativeFilenames::combinePathSegments(turbomoleExecutableBase_, binaryName);

  if (outputToFile) {
    if (outputFileName.empty()) {
      outputFileName = binaryName + ".out";
    }
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
  bp::ipstream _stderr;
  std::string binary = NativeFilenames::combinePathSegments(turbomoleExecutableBase_, binaryName);
  bp::child c(binary, bp::std_in<stdInFile, bp::std_out> bp::null, bp::std_err > _stderr, workingDir);
  c.wait(); // Wait for the process to exit

  bool success = jobWasSuccessful(_stderr);
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
  auto it = std::find(std::begin(basisSetNotationsLower_), std::end(basisSetNotationsLower_), basisSetString);
  if (it != basisSetNotationsLower_.end()) {
    assert(std::distance(basisSetNotationsLower_.begin(), it) < static_cast<long>(basisSetNotations_.size()));
    basisSetString = basisSetNotations_[std::distance(basisSetNotationsLower_.begin(), it)];
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
