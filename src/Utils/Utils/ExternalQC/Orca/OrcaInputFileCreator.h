/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_EXTERNALQC_ORCAINPUTFILECREATOR_H
#define UTILS_EXTERNALQC_ORCAINPUTFILECREATOR_H

#include <ostream>
#include <string>

namespace Scine {
namespace Utils {

class AtomCollection;
class PropertyList;
class Settings;

namespace ExternalQC {

/**
 * @class OrcaInputFileCreator OrcaInputFileCreator.h
 * @brief This class creates ORCA input files.
 */
class OrcaInputFileCreator {
 public:
  /**
   * @brief Create the ORCA input file with the filename 'filename'
   * @param filename Name of the input file.
   * @param atoms Molecular structure.
   * @param settings Settings of the calculation.
   * @param requiredProperties All the required properties (energy, gradients,...)
   */
  void createInputFile(const std::string& filename, const AtomCollection& atoms, const Settings& settings,
                       const PropertyList& requiredProperties);
  /**
   * @brief Create the ORCA input file with the filename 'filename'
   * @param out std::ostream.
   * @param atoms Molecular structure.
   * @param settings Settings of the calculation.
   * @param requiredProperties All the required properties (energy, gradients,...)
   */
  void createInputFile(std::ostream& out, const AtomCollection& atoms, const Settings& settings,
                       const PropertyList& requiredProperties);

 private:
  void printCalculationType(std::ostream& out, const Settings& settings, const PropertyList& requiredProperties);
  void printTitle(std::ostream& out);
  void printStructure(std::ostream& out, const AtomCollection& atoms, const Settings& settings);
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_EXTERNALQC_ORCAINPUTFILECREATOR_H