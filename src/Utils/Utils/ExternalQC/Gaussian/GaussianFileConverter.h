/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_GAUSSIANFILECONVERTER_H
#define UTILS_GAUSSIANFILECONVERTER_H

namespace Scine {
namespace Utils {
namespace ExternalQC {

/**
 * @brief Functionality to convert Gaussian checkpoint files to formatted Gaussian checkpoint files and vice versa
 *
 */
namespace GaussianFileConverter {

/**
 * @brief Generate a formatted checkpoint file from a binary checkpoint file
 *
 * @param fileBase The base of the name of the chk file and the new fchk file
 * @param workingDirectory The directory with the files
 * @param gaussianDirectory The directory the "formchk" executable is located in.
 */
std::string generateFormattedCheckpointFile(const std::string& fileBase, const std::string& workingDirectory,
                                            const std::string& gaussianDirectory);

/**
 * @brief Generate a binary checkpoint file from a formatted checkpoint file
 *
 * @param fileBase The base of the name of the fchk file and the new chk file
 * @param workingDirectory The directory with the files
 * @param gaussianDirectory The directory the "unchk" executable is located in.
 */
std::string generateCheckpointFile(const std::string& fileBase, const std::string& workingDirectory,
                                   const std::string& gaussianDirectory);

/**
 * @brief A custom exception for non-existing checkpoint files or formatted checkpoint files
 */
class GaussianCheckpointFileNotFoundException : public std::runtime_error {
 public:
  GaussianCheckpointFileNotFoundException(const std::string& message) : std::runtime_error(message.c_str()) {
  }
  const char* what() const noexcept override {
    return std::runtime_error::what();
  }
};

class GaussianCheckpointConversionException : public std::runtime_error {
 public:
  GaussianCheckpointConversionException(const std::string& message) : std::runtime_error(message.c_str()) {
  }
  const char* what() const noexcept override {
    return std::runtime_error::what();
  }
};

} // namespace GaussianFileConverter
} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_GAUSSIANFILECONVERTER_H
