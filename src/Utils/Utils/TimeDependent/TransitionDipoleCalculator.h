/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_TRANSITIONDIPOLECALCULATOR_H
#define UTILS_TRANSITIONDIPOLECALCULATOR_H

#include <Utils/DataStructures/DipoleMatrix.h>
#include <Utils/Math/IterativeDiagonalizer/SpinAdaptedEigenContainer.h>
#include <Utils/Scf/LcaoUtils/ElectronicOccupation.h>
#include <Utils/Typenames.h>
#include <exception>

namespace Scine {
namespace Utils {
namespace LcaoUtils {

class ElectronicOccupation;
} // namespace LcaoUtils

class DipoleMatrix;

struct Excitation {
  int occ;
  int vir;
  auto operator==(const Excitation& lhs) const -> bool {
    return (occ == lhs.occ) && (vir == lhs.vir);
  }
  auto operator!=(const Excitation& lhs) const -> bool {
    return (occ != lhs.occ) || (vir != lhs.vir);
  }
};
/**
 * @brief Class to compute the transition dipole moment and the corresponding oscillator strength.
 * The transition dipole moment and the corresponding oscillator strength
 * in molecular orbital basis are calculated according to
 * J. D. Baker, M. C. Zerner,
 * Applications of the random phase approximation with the INDO/S Hamiltonian: UV-VIS spectra of free base porphin,
 * Chem. Phys. Lett., 1990, 175, 192-196
 * as
 * \f$ d_v = \sqrt{2} * \sum_{i,j} D_{ij} * X_{ij,v} \f$
 * and
 * \f$ s_v = 2/3 * E_v^{-1} * d_v^2 \f$
 * in restricted formalism. For an unrestricted reference, the transition dipole moment is calculated as
 * \f$ d_v = \sum_{i,j} D_{ij} * X_{ij,v} \f$
 *
 * d_v: transition dipole vector.
 * s_v: oscillator strength.
 * E_v: transition energy of the v-th state.
 * i,j: molecular orbitals.
 * X_v: transition state vector.
 *
 */
struct TransitionDipoleCalculator {
  /**
   * @brief Calculates the x component of the transition dipole moment.
   * @param dipoleMatrix A dipole moment matrix containing the dipole integrals. In MO basis.
   * @param stateVector A state vector as recovered from a CIS/TD-DFT calculation.
   * @param excitations Vector containing the hole-particle excitations.
   * @return The x component of the transition dipole vector.
   */
  template<Utils::Reference restrictedness>
  static double x(const Utils::DipoleMatrix& dipoleMatrix, const Eigen::VectorXd& stateVector,
                  const std::vector<Excitation>& excitations) {
    return evaluate<restrictedness>(dipoleMatrix.x().get<Utils::DerivativeOrder::Zero>(), stateVector, excitations);
  }
  /**
   * @brief Calculates the x component of the transition dipole moment.
   * @param dipoleMatrix A dipole moment matrix containing the dipole integrals. In MO basis.
   * @param stateVector A state vector as recovered from a CIS/TD-DFT calculation.
   * @param excitations Vector containing the hole-particle excitations.
   * @return The x component of the transition dipole vector.
   */
  template<Utils::Reference restrictedness>
  static double y(const Utils::DipoleMatrix& dipoleMatrix, const Eigen::VectorXd& stateVector,
                  const std::vector<Excitation>& excitations) {
    return evaluate<restrictedness>(dipoleMatrix.y().get<Utils::DerivativeOrder::Zero>(), stateVector, excitations);
  }
  /**
   * @brief Calculates the x component of the transition dipole moment.
   * @param dipoleMatrix A dipole moment matrix containing the dipole integrals. In MO basis.
   * @param stateVector A state vector as recovered from a CIS/TD-DFT calculation.
   * @param excitations Vector containing the hole-particle excitations.
   * @return The x component of the transition dipole vector.
   */
  template<Utils::Reference restrictedness>
  static double z(const Utils::DipoleMatrix& dipoleMatrix, const Eigen::VectorXd& stateVector,
                  const std::vector<Excitation>& excitations) {
    return evaluate<restrictedness>(dipoleMatrix.z().get<Utils::DerivativeOrder::Zero>(), stateVector, excitations);
  }
  /**
   * @brief Calculates the transition dipole moment vector.
   * TODO: If a triplet state is calculated, this does not give zero as the state
   * is just multiplied by 2 (as should be in singlet state), but nothing is done for the triplet
   * state. (Combination in singlet ->
   * |CSF^1> = 1/sqrt(2)(|alpha> + |beta>) -> 2 * spatial part
   * |CSF^3> = 1/sqrt(2)(|alpha> - |beta>) -> 0 * spatial part
   * @param dipoleMatrix A dipole moment matrix containing the dipole integrals. In MO basis.
   * @param stateVectors The state vectors as recovered from a CIS/TD-DFT calculation.
   * @param excitations Vector containing the hole-particle excitations.
   * @return The transition dipole vector for all the states.
   */
  template<Utils::Reference restrictedness>
  static Eigen::Matrix3Xd calculate(const Utils::DipoleMatrix& dipoleMatrix, const Eigen::MatrixXd& stateVectors,
                                    const std::vector<Excitation>& excitations) {
    Eigen::Matrix3Xd transitionDipoles(3, stateVectors.cols());
    for (int i = 0; i < stateVectors.cols(); ++i) {
      Utils::Dipole transitionDipole(3);
      transitionDipole << x<restrictedness>(dipoleMatrix, stateVectors.col(i), excitations),
          y<restrictedness>(dipoleMatrix, stateVectors.col(i), excitations),
          z<restrictedness>(dipoleMatrix, stateVectors.col(i), excitations);
      transitionDipoles.col(i) = transitionDipole;
    }
    return transitionDipoles;
  }
  /**
   * @brief Transforms the transition dipole moment vector to the corresponding oscillator strength.
   * @param transitionDipoleMoments The transition dipole moments.
   * @param eigenvalues The transition energies for the calculated transitions.
   * @return A Eigen::VectorXd containing a single oscillator strength for each state.
   */
  static Eigen::VectorXd transitionDipoleMomentToOscillatorStrength(const Eigen::Matrix3Xd& transitionDipoleMoments,
                                                                    const Eigen::VectorXd& eigenvalues) {
    assert(transitionDipoleMoments.cols() == eigenvalues.size());
    Eigen::VectorXd transitionDipoleSquared = transitionDipoleMoments.cwiseProduct(transitionDipoleMoments).colwise().sum();
    return 2.0 / 3.0 * transitionDipoleSquared.cwiseProduct(eigenvalues);
  }

 private:
  template<Utils::Reference restrictedness>
  static double evaluate(const Eigen::MatrixXd& dipoleComponentMatrix, const Eigen::VectorXd& stateVector,
                         const std::vector<Excitation>& excitations) {
    int nExc = excitations.size();
    assert(nExc == stateVector.size());
    double sum = 0.0;
    for (int index = 0; index < nExc; ++index) {
      int occ = excitations[index].occ;
      int vir = excitations[index].vir;
      sum += dipoleComponentMatrix(occ, vir) * stateVector[index];
    }
    return factor<restrictedness>() * sum;
  }
  //! @brief 1.0 for unrestricted, sqrt(2) for restricted.
  template<Utils::Reference restrictedness>
  static double factor();
};

template<>
inline auto TransitionDipoleCalculator::factor<Utils::Reference::Restricted>() -> double {
  return std::sqrt(2);
}
template<>
inline auto TransitionDipoleCalculator::factor<Utils::Reference::Unrestricted>() -> double {
  return 1.0;
}
} // namespace Utils
} // namespace Scine

#endif // UTILS_TRANSITIONDIPOLECALCULATOR_H
