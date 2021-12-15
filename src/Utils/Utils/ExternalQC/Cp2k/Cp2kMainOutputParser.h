/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_EXTERNALQC_CP2KMAINOUTPUTPARSER_H
#define UTILS_EXTERNALQC_CP2KMAINOUTPUTPARSER_H

#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/DataStructures/AtomsOrbitalsIndexes.h>
#include <Utils/DataStructures/DensityMatrix.h>
#include <Utils/Scf/LcaoUtils/SpinMode.h>
#include <Utils/Typenames.h>
#include <Eigen/Core>
#include <regex>
#include <string>

namespace Scine {
namespace Utils {
namespace ExternalQC {

/**
 * @class Cp2kMainOutputParser Cp2kMainOutputParser.h
 * @brief This class parses information out of the main CP2K output file.
 */
class Cp2kMainOutputParser {
 public:
  /**
   * @brief Constructor.
   * @param outputFileName Name of the output file.
   */
  explicit Cp2kMainOutputParser(const std::string& outputFileName, const std::string& additionalOutput = "");
  /**
   * @brief Parse the CP2K output for errors, it is expected to do nothing if no error is present
   * @throws OutputFileParsingError if error is present
   */
  void checkForErrors() const;
  /**
   * @brief Parse the energy from the CP2K output.
   * @return The energy as a double.
   */
  double getEnergy() const;
  /**
   * @brief Parse the gradients from the CP2K output.
   * @return GradientCollection
   */
  GradientCollection getGradients() const;
  /**
   * @brief Parse the Hessian matrix from the output file.
   * @return The Hessian matrix.
   */
  HessianMatrix getHessian() const;
  /**
   * @brief Parse the Hirshfeld charges from the CP2K output.
   * @return The Hirshfeld charges.
   */
  std::vector<double> getHirshfeldCharges() const;
  /**
   * @brief Parse the number of gaussian functions per grid from the CP2K output.
   * @return The numbers for each grid starting with the highest energy cutoff grid.
   */
  std::vector<int> getGridCounts() const;
  /**
   * @brief Parse different matrices from the CP2K output and construct Mayer bond orders from that.
   * @param The elements of the structure.
   * @param The spin mode such as restricted or unrestricted.
   * @return The Mayer bond orders.
   */
  BondOrderCollection getBondOrders(const ElementTypeCollection& elements, SpinMode spinMode) const;
  /**
   * @brief Parse the density matrix from the CP2K output.
   * @param The spin mode such as restricted or unrestricted.
   * @return The density matrix.
   */
  DensityMatrix getDensityMatrix(SpinMode spinMode) const;
  /**
   * @brief Parse the overlap matrix from the CP2K output.
   * @return The overlap matrix.
   */
  Eigen::MatrixXd getOverlapMatrix() const;
  /**
   * @brief Parse the number of spherical basis functions for each element from the CP2K output and construct and Atom
   * to AO index from this information.
   * @param The elements of the structure.
   * @return The index specifying the connection between atoms and AOs.
   */
  AtomsOrbitalsIndexes getAtomAoIndex(const ElementTypeCollection& elements) const;
  /**
   * @brief Parse the number of total AOs from the CP2K output.
   * @return The number of AOs.
   */
  int getNumberOfAos() const;
  /**
   * @brief Parse the number of electrons from the CP2K output.
   * @note in case of a restricted calculation the result contains only one entry, the total number of electrons; in
   * case of an unrestricted calculation, the number of alpha and beta electrons are stored as the first and second
   * entry.
   * @return The number of electrons.
   */
  std::vector<int> getNumberOfElectrons() const;
  /**
   * @brief Parse the molecular symmetry for which the thermochemistry should be computed.
   * @return The molecular symmetry number.
   */
  int getSymmetryNumber() const;

 private:
  // @brief extract the last matrix print out block from a CP2K output file content
  std::string parseMatrixStringBlockFromFileContent(const std::string& content, const std::regex& startKeyword) const;
  // @brief extract a matrix from a type of CP2K matrix print out block
  Eigen::MatrixXd parseMatrixFromStringBlock(const std::string& block, const std::string& propName, int dimension) const;
  static std::string extractContent(const std::string& filename);
  void extractRuntype();

  std::string content_;
  std::string additionalContent_;
  std::string runtype_;
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_EXTERNALQC_CP2KMAINOUTPUTPARSER_H