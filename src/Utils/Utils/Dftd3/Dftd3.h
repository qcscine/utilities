/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SCINE_DFTD3_H
#define SCINE_DFTD3_H

#include "Dftd3Atom.h"
#include "Dftd3Parameters.h"
#include "Utils/Math/AtomicSecondDerivativeCollection.h"
#include "Utils/Math/DerivOrderEnum.h"
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Math/AutomaticDifferentiation/AutomaticDifferentiationHelpers.h>
#include <Utils/Math/AutomaticDifferentiation/Second1D.h>
#include <Utils/Typenames.h>
#include <vector>

namespace Scine {
namespace Utils {
namespace Tests {
class ADftd3Calculation_EnergyCalculationIsInternallyConsistent_Test;
class ADftd3Calculation_HasCorrectDerivativeOfC6WrtCoordNumber_Test;
class ADftd3Calculation_HasCorrectDerivativeOfEnergyWrtDistance_Test;
} // namespace Tests
namespace Dftd3 {
/**
 * @class Dftd3 Dftd3.h
 * @brief A class calculating D3 semi-classical dispersion corrections (with BJ damping) for energies and gradients.
 *
 *      This class uses a default constructor. The calculations are initialized by the initialize function.
 *      Simple calculations are done using the calculate function, to which you can pass the required
 *      derivative level.
 *
 * References:
 *      - J. Chem. Phys. 132, 154104, (2010)  (Description D3)
 *      - J. Comput. Chem. 32, 1456, (2011)   (Description of BJ damping function)
 *
 * Developer Info:
 *      - The assignment of some functions as public is based on the fact that they are needed
 *        for the non-bonded interactions in the molecular mechanics library.
 */
class Dftd3 {
 public:
  /**
   * @brief The initialize function needs to be called before a calculation as a set-up.
   *
   *        The meaning of the three Grimme D3 parameters can be found in the aformentioned reference.
   */
  void initialize(AtomCollection atomCollection, double a1, double s8, double a2);
  /**
   * @brief Getter for the D3 energy.
   * @return double
   */
  double getEnergy();
  /**
   * @brief Getter for the D3 gradients.
   */
  GradientCollection getGradients();
  /**
   * @brief Getter for the structure in the format of a vector of Dftd3Atoms.
   */
  std::vector<Dftd3Atom> getStructure();
  /**
   * @brief This function needs to be called to do a D3 dispersion correction calculation.
   * @param d derivative type, up to which will be calculated (either none or first).
   */
  void calculate(derivativeType d);
  /**
   * @brief Set structure as vector of Dftd3Atoms.
   */
  void setStructure(std::vector<Dftd3Atom> structure);
  /**
   * @brief Calculates C6 coefficient for one atom pair.
   * @return double
   */
  double calculateC6Coefficient(Dftd3Atom atom1, Dftd3Atom atom2);
  /**
   * @brief Calculates C8 coefficient for one atom pair.
   * @param C6 coefficient for that atom pair.
   * @return double
   */
  double calculateC8Coefficient(Dftd3Atom atom1, Dftd3Atom atom2, double c6);
  /**
   * @brief Getter for the pair cutoff radius R0 for one atom pair.
   * @return double
   */
  double getR0(int atom1Index, int atom2Index);
  /**
   * @brief Getter for the matrix of all C6 coefficients.
   * @return Eigen::MatrixXd
   */
  Eigen::MatrixXd getC6Matrix() const;
  /**
   * @brief This function calculates the D3 coordination number (fractional value).
   * @return double
   */
  double calculateCoordinationNumber(Dftd3Atom atom);

 private:
  /*
   * @brief Declare some googletests as friends in order to make it possible to test private functions.
   */
  friend class Tests::ADftd3Calculation_EnergyCalculationIsInternallyConsistent_Test;
  friend class Tests::ADftd3Calculation_HasCorrectDerivativeOfC6WrtCoordNumber_Test;
  friend class Tests::ADftd3Calculation_HasCorrectDerivativeOfEnergyWrtDistance_Test;
  /*
   * @brief Energy evaluator for one atom pair.
   * @return double
   */
  double evaluateEnergy(Dftd3Atom& atom1, Dftd3Atom& atom2);
  /*
   * @brief Gradient evaluator for one atom pair.
   * @param dE_drij This is the simple derivative of the energy with respect to the distance rij.
   * @param dc6 A more complicated vector of doubles need for the computation of derivatives.
   *        Definition is given in the implementation (.cpp).
   * @param UpToSecondDerivatives The container holds the current derivative. This container can in principle hold
   *        all values up to the second derivatives. It is given by reference and modified within this function.
   *        Thus, there is no return value.
   */
  void evaluateGradients(Dftd3Atom& atom1, Dftd3Atom& atom2, double& dE_drij, std::vector<double>& dc6,
                         AtomicSecondDerivativeCollection& UpToSecondDerivatives);
  /*
   * @brief This function calculates the matrices that hold all C6 and C8 coefficients as well as R0 cutoff radii.
   */
  void calculateValuesForConstants();
  /*
   * @brief Evaluation of the derivative of the energy of one atom pair w.r.t the C6 coefficient.
   *
   *        The derivative dE/dC6 ist just equal to E/C6, therefore the name of the function.
   *
   * @return double
   */
  double evaluateEnergyDividedByC6(Dftd3Atom& atom1, Dftd3Atom& atom2);
  /*
   * @brief Evaluation of the derivative of the coord. number of atom 1 or atom 2 (that is the same) w.r.t. the distance
   *        between atom 1 and atom 2.
   * @return double
   */
  double evaluateGradientOfCoordNumberWrtDistance(Dftd3Atom& atom1, Dftd3Atom& atom2);
  /*
   * @brief Evaluation of the energy gradients with respect to the distance between atom 1 and atom 2.
   * @return double derivative d(E)/d(r12)
   */
  double evaluateGradientsWrtDistances(Dftd3Atom& atom1, Dftd3Atom& atom2);
  /*
   * @brief Evaluation of the derivative of the C6 coefficient of atom 1 and atom 2 w.r.t. the coord. number of atom 1.
   * @return double
   */
  double evaluateGradientOfC6WrtCoordNumber(Dftd3Atom& atom1, Dftd3Atom& atom2);

  std::vector<Dftd3Atom> structure_;
  double energy_;
  GradientCollection gradients_;
  // @brief container for derivatives.
  AtomicSecondDerivativeCollection UpToSecondDerivatives_;
  // @brief D3 parameters defined in Dftd3Parameters.h
  Dftd3Parameters parameters_;
  // @brief Matrix holding the derivatives of the C6 coefficients w.r.t a coordination number.
  Eigen::MatrixXd dC6_dCN_;
  // @brief Matrix holding the derivatives of the coord. numbers w.r.t pair distances rij.
  Eigen::MatrixXd dCN_drij_;
  // @brief All C6 coefficients.
  Eigen::MatrixXd C6_;
  // @brief All C8 coefficients.
  Eigen::MatrixXd C8_;
  // @brief All R0 cutoff radii.
  Eigen::MatrixXd R0_;
};

} // namespace Dftd3
} // namespace Utils
} // namespace Scine

#endif // SCINE_DFTD3_H
