/**
 * @file
 * @brief Default module that is always loaded by Core
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILSOS_MRCCHELPER_H
#define UTILSOS_MRCCHELPER_H

#include <string>

namespace Scine::Utils::ExternalQC {
/**
 * @class
 * @brief A struct to organize the files used with MRCC.
 */
struct MrccFiles {
  /**
   * @brief Constructor
   * @param baseDir The directory in which the calculations should be run.
   */
  MrccFiles(const std::string& baseDir);

  std::string input;
  std::string output;
};

/**
 * @class
 * @brief Helper to run MRCC calculations.
 */
class MrccHelper {
 public:
  MrccHelper(const std::string& mrccBasePath, const std::string& workingDirectory);

  const MrccFiles& getFiles();
  /**
   * @brief Run MRCC as `dmrcc > dmrcc.out`.
   */
  void run();

 private:
  const std::string mrccBasePath_;
  const std::string workingDirectory_;
  const std::string dmrccExecutable_;
  const std::string ccsdExecutable_;
  const std::string scfExecutable_;
  const MrccFiles files_;
};

namespace MrccExecutableNames {
static const std::string dmrcc = "dmrcc";
static const std::string ccsd = "ccsd";
static const std::string scf = "scf";
} // namespace MrccExecutableNames

namespace MrccFileNames {
/// @brief The input file name is fixed by MRCC to MINP.
static const std::string input = "MINP";
/// @brief Default output file name for MRCC.
static const std::string output = "dmrcc.out";
} // namespace MrccFileNames

} // namespace Scine::Utils::ExternalQC

#endif // UTILSOS_MRCCHELPER_H
