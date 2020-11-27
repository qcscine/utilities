/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_EXTERNALQC_ORCAMAINOUTPUTPARSER_H
#define UTILS_EXTERNALQC_ORCAMAINOUTPUTPARSER_H

#include <Utils/Typenames.h>
#include <Eigen/Core>
#include <string>

namespace Scine {
namespace Utils {
class BondOrderCollection;
namespace ExternalQC {

/**
 * @class OrcaMainOutputParser OrcaMainOutputParser.h
 * @brief This class parses information out of the main ORCA output file.
 */
class OrcaMainOutputParser {
 public:
  /**
   * @brief Constructor.
   * @param outputFileName Name of the output file.
   */
  explicit OrcaMainOutputParser(const std::string& outputFileName);
  /**
   * @brief Parse the ORCA output for errors, it is expected to do nothing if no error is present
   * @throws OutputFileParsingError if error is present
   */
  void checkForErrors() const;
  /**
   * @brief Parse the energy from the ORCA output.
   * @return The energy as a double.
   */
  double getEnergy() const;
  /**
   * @brief Parse the Mayer bond orders from the ORCA output.
   * @return The Utils::BondOrderCollection.
   */
  Utils::BondOrderCollection getBondOrders() const;
  /**
   * @brief Parse the Hirshfeld charges from the ORCA output.
   * @return The Hirshfeld charges.
   */
  std::vector<double> getHirshfeldCharges() const;
  /**
   * @brief Parse the number of atoms from the ORCA output.
   * @return The number of atoms.
   */
  int getNumberAtoms() const;
  /**
   * @brief Parse the gradients from the ORCA output.
   * @return GradientCollection
   */
  GradientCollection getGradients() const;
  /**
   * @brief Parse the temperature for which the thermochemistry was computed.
   * @return The temperature in Kelvin.
   */
  double getTemperature() const;
  /**
   * @brief Parse the enthalpy from the output.
   * @return The enthalpy in Hartree.
   */
  double getEnthalpy() const;
  /**
   * @brief Parse the entropy from the output.
   * @return The entropy in Hartree/Kelvin.
   */
  double getEntropy() const;
  /**
   * @brief Parse the zero point vibrational energy from the output.
   * @return The zero point vibrational energy in Hartree.
   */
  double getZeroPointVibrationalEnergy() const;
  /**
   * @brief Parse the Gibbs free energy from the output.
   * @return The Gibbs free energy in Hartree.
   */
  double getGibbsFreeEnergy() const;

  /**
   * @brief Parse the molecular symmetry for which the thermochemistry was computed.
   * @return The molecular symmetry number.
   */
  double getSymmetryNumber() const;

 private:
  void extractContent(const std::string& filename);

  std::string content_;
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_EXTERNALQC_ORCAMAINOUTPUTPARSER_H