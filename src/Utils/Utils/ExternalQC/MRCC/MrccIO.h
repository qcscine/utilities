/**
 * @file
 * @brief Default module that is always loaded by Core
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILSOS_MRCCIO_H
#define UTILSOS_MRCCIO_H

#include <Core/BaseClasses/ObjectWithLog.h>
#include <Utils/ExternalQC/MRCC/MrccHelper.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Settings.h>
#include <fstream>
#include <string>

namespace Scine {
namespace Utils {
namespace ExternalQC {
/**
 * @class
 * @brief The supported MRCC methods.
 */
enum class MrccMethod { HF = 0, DFT = 1, LNO_MP2 = 2, LNO_CCSD = 3, LNO_CCSD_T = 4 };
/**
 * @class
 * @brief This class handles the input file generation and output file parsing for MRCC.
 */
class MrccIO : public Core::ObjectWithLog {
 public:
  /**
   * @brief Constructor.
   * @param files         The MRCC files, e.g., input and output.
   * @param settings      The calculation settings.
   * @param methodFamily  The method family, e.g., HF, DFT, MP2, or CC.
   */
  MrccIO(const MrccFiles& files, const Settings& settings, const std::string& methodFamily);

  /**
   * @brief Write the MRCC input file for the given molecule and settings.
   * @param atoms The molecular structure.
   */
  void writeInput(const AtomCollection& atoms);
  /**
   * @brief Read the output file.
   */
  std::string readOutput();
  /**
   * @brief Getter for the energy.
   * @param output The output file as a string.
   * @return The energy.
   */
  double getEnergy(const std::string& output);

 private:
  const MrccFiles files_;
  const Settings settings_;
  const MrccMethod mrccMethod_;

  ///@brief Resolve the method to the possible MrccMethod-options.
  static MrccMethod getMrccMethod(const Settings& settings, const std::string& methodFamily);
  void ensureSuccessFullCalculation(const std::string& out);

  // Add keyword blocks to the MRCC input
  void addCoordinateDefinition(const AtomCollection& atoms, std::ofstream& outStream);
  void addAllowedResources(std::ofstream& outStream);
  void addChargeAndMultiplicity(std::ofstream& outStream);
  void addMethodDefinition(std::ofstream& outStream);
  void addSCFKeywords(std::ofstream& outStream);
  void addCalcKeyword(std::ofstream& outStream);
  void addSolvationKeywords(std::ofstream& outStream);
  void addBasisSetKeyword(std::ofstream& outStream);
  void addSCFTypeKeyword(std::ofstream& outStream);
  void addLocalCorrelationKeywords(std::ofstream& outStream);

  bool isLocalCorrelation();
  std::string getLNOThresholdsFromMethod();
  std::string getEnergyString();
  std::string functionalInMrccFormat();
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILSOS_MRCCIO_H
