/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_EXTERNALQC_EXTERNALPROGRAM_H
#define UTILS_EXTERNALQC_EXTERNALPROGRAM_H

#include <Utils/Geometry/AtomCollection.h>
#include <string>

namespace Scine {
namespace Utils {
namespace ExternalQC {

/**
 * @class ExternalProgram ExternalProgram.h
 * @brief This class allows for running external programs through SCINE.
 *
 *        This class is used by the ExternalQC calculators.
 */
class ExternalProgram {
 public:
  /**
   * @brief Setter for the working directory.
   */
  void setWorkingDirectory(const std::string& directory);
  /**
   * @brief Getter for the working directory.
   */
  const std::string& getWorkingDirectory() const;

  /**
   * @brief Create a temporary working directory.
   */
  void createWorkingDirectory();
  /**
   * @brief Execute a command with no input or output file.
   */
  void executeCommand(const std::string& command) const;
  /**
   * @brief Execute a command with an output file given as an argument.
   */
  void executeCommand(const std::string& command, const std::string& outputFile) const;
  /**
   * @brief Execute a comand with an input and output file given as arguments.
   */
  void executeCommand(const std::string& command, const std::string& inputFile, const std::string& outputFile) const;
  /**
   * @brief Generate filename by prepending the working directory.
   */
  std::string generateFullFilename(const std::string& filename) const;
  /**
   * @brief Getter of file for standard error output.
   */
  std::string getErrorOutFile() const;
  /**
   * @brief Setter to write error out to specific file.
   */
  void setErrorOutFile(const std::string& filename);

 private:
  int executeCommandImpl(const std::string& command, const std::string& inputFile, const std::string& outputFile) const;
  std::string workingDirectory_{};
  std::string _errorOut;
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_EXTERNALQC_EXTERNALPROGRAM_H
