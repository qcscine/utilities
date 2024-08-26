/**
 * @file
 * @brief Default module that is always loaded by Core
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "MrccHelper.h"
#include "Utils/IO/NativeFilenames.h"
#include <boost/process.hpp>
#include <regex>

namespace Scine {
namespace Utils {
namespace ExternalQC {

MrccHelper::MrccHelper(const std::string& mrccBasePath, const std::string& workingDirectory)
  : mrccBasePath_(mrccBasePath),
    workingDirectory_(workingDirectory),
    dmrccExecutable_(NativeFilenames::combinePathSegments(mrccBasePath_, MrccExecutableNames::dmrcc)),
    ccsdExecutable_(NativeFilenames::combinePathSegments(mrccBasePath_, MrccExecutableNames::ccsd)),
    scfExecutable_(NativeFilenames::combinePathSegments(mrccBasePath_, MrccExecutableNames::scf)),
    files_(workingDirectory_) {
  if (!boost::filesystem::exists(this->dmrccExecutable_)) {
    throw std::runtime_error("MRCC calculator cannot locate the dmrcc executable at " + this->dmrccExecutable_);
  }
  if (!boost::filesystem::exists(this->ccsdExecutable_)) {
    throw std::runtime_error("MRCC calculator cannot locate the ccsd executable at " + this->ccsdExecutable_);
  }
  if (!boost::filesystem::exists(this->scfExecutable_)) {
    throw std::runtime_error("MRCC calculator cannot locate the scf executable at " + this->scfExecutable_);
  }
}

const MrccFiles& MrccHelper::getFiles() {
  return this->files_;
}

void MrccHelper::run() {
  const auto workingDirectory = boost::process::start_dir(this->workingDirectory_);
  boost::process::ipstream _stderr;
  boost::filesystem::remove(files_.output);
  boost::process::child c(this->dmrccExecutable_, boost::process::std_out > files_.output,
                          boost::process::std_err > _stderr, workingDirectory);
  c.wait();
}

MrccFiles::MrccFiles(const std::string& baseDir) {
  input = NativeFilenames::combinePathSegments(baseDir, MrccFileNames::input);
  output = NativeFilenames::combinePathSegments(baseDir, MrccFileNames::output);
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
