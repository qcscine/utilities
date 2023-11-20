/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_EXTERNALQC_GAUSSIANINPUTFILECREATOR_H
#define UTILS_EXTERNALQC_GAUSSIANINPUTFILECREATOR_H

#include <ostream>
#include <string>

namespace Scine {
namespace Utils {

class AtomCollection;
class PropertyList;
class Settings;

namespace ExternalQC {

/**
 * @class GaussianInputFileCreator GaussianInputFileCreator.h
 * @brief This class creates Gaussian input files.
 */
class GaussianInputFileCreator {
 public:
  /**
   * @brief Create the Gaussian input file with the filename 'filename'
   * @param filename Name of the input file.
   * @param checkpointFilename Name of the checkpoint file.
   * @param atoms Molecular structure.
   * @param settings Settings of the calculation.
   * @param requiredProperties All the required properties (energy, gradients,...)
   */
  static void createInputFile(const std::string& filename, const std::string& checkpointFilename, const AtomCollection& atoms,
                              const Settings& settings, const PropertyList& requiredProperties);
  /**
   * @brief Create the Gaussian input file with the filename 'filename'
   * @param out std::ostream.
   * @param checkpointFilename Name of the checkpoint file.
   * @param atoms Molecular structure.
   * @param settings Settings of the calculation.
   * @param requiredProperties All the required properties (energy, gradients,...)
   */
  static void createInputFile(std::ostream& out, const std::string& checkpointFilename, const AtomCollection& atoms,
                              const Settings& settings, const PropertyList& requiredProperties);

 private:
  static void printCalculationType(std::ostream& out, const std::string& checkpointFilename, const Settings& settings,
                                   const PropertyList& requiredProperties);
  static void printTitle(std::ostream& out);
  static void printStructure(std::ostream& out, const AtomCollection& atoms, const Settings& settings);
  static inline std::string translateDispersion(const std::string& scineDispersion) {
    return (scineDispersion.empty()) ? "" : "EmpiricalDispersion=G" + scineDispersion;
  };
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_EXTERNALQC_GAUSSIANINPUTFILECREATOR_H
