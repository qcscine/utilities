/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_CONCEPTUALDFT_H
#define UTILS_CONCEPTUALDFT_H

#include <Eigen/Core>

namespace Scine {
namespace Utils {

/**
 * @brief Functionality to calculate a set of conceptual DFT properties based on the finite difference
 * approximation.
 *
 * To use this functionality the energy and/or atomic charges of the optimized structure of interest and of the same
 * structure (not reoptimized!) with one additional electron and one electron less are needed.
 *
 * An introduction into conceptual DFT can be found for example in:
 * Chattaraj, P. K. Chemical Reactivity Theory : A Density Functional View; CRC Press, 2009.
 * https://doi.org/10.1201/9781420065442.
 */
namespace ConceptualDft {
/**
 * @brief Struct containing global conceptual DFT parameters.
 */
struct GlobalConceptualDftContainer {
  double chemicalPotential; /**<\f$\mu=\frac{E(R, N+1) - E(R, N-1)}{2}\f$*/
  double electronegativity; /**<\f$\chi=\frac{E(R, N-1) - E(R, N+1)}{2}\f$*/
  double hardness;          /**<\f$\eta=E(R, N-1)+E(R, N+1)-2E(R, N)\f$*/
  double softness;          /**<\f$S=\frac{1}{\eta}\f$*/
  double electrophilicity;  /**<\f$\omega=\frac{\mu^2}{2\eta}\f$*/
};

/**
 * @brief Struct containing local conceptual DFT properties.
 */
struct LocalConceptualDftContainer {
  Eigen::VectorXd fukuiPlus;      /**<\f$f^+_k=q_k(R, N) - q_k(R, N+1)\f$*/
  Eigen::VectorXd fukuiMinus;     /**<\f$f^-_k=q_k(R, N-1) - q_k(R, N)\f$*/
  Eigen::VectorXd fukuiRadical;   /**<\f$f^0_k=\frac{1}{2}(q_k(R, N-1) - q_k(R, N+1))\f$*/
  Eigen::VectorXd dualDescriptor; /**<\f$\Delta f_k = f_k^+ - f_k^-\f$*/
};

/**
 * @brief Struct containing local and global conceptual DFT properties.
 */
struct ConceptualDftContainer {
  GlobalConceptualDftContainer global; /**< The global indices */
  LocalConceptualDftContainer local;   /**< The condensed to atom indices */
};
/**
 * @brief Calculates a set of global and local conceptual DFT parameters.
 *
 * @note The quality of the resulting Fukui and dual descriptor indices heavily depends on the
 * quality of the atomic charges. We recommend Hirshfeld charges.
 *
 * @param energy The energy of the structure of interest \f$E(R, N)\f$.
 * @param atomicCharges The atomic charges of the structure of interest \f$q(R, N)\f$.
 * @param energyPlus The energy for the same geometry with one additional electron \f$E(R, N+1)\f$
 * @param atomicChargesPlus The atomic charges for the same geometry with one additional electron \f$q(R, N+1)\f$.
 * @param energyMinus The energy for the same geometry with one electron less \f$E(R, N-1)\f$
 * @param atomicChargesMinus The atomic charges for the same geometry with one electron less \f$q(R, N-1)\f$.
 * @return A ConceptualDftContainer with local and global cDFT properties.
 */
ConceptualDftContainer calculate(const double energy, const Eigen::VectorXd& atomicCharges, const double energyPlus,
                                 const Eigen::VectorXd& atomicChargesPlus, const double energyMinus,
                                 const Eigen::VectorXd& atomicChargesMinus);

/**
 * @brief Calculates a set of global conceptual DFT parameters.
 *
 * @param energy The energy of the structure of interest \f$E(R, N)\f$.
 * @param energyPlus The energy for the same geometry with one additional electron \f$E(R, N+1)\f$
 * @param energyMinus The energy for the same geometry with one electron less \f$E(R, N-1)\f$
 * @return A GlobalConceptualDftContainer with global cDFT properties.
 */
GlobalConceptualDftContainer calculateGlobal(const double energy, const double energyPlus, const double energyMinus);

/**
 * @brief Calculates the condensed to atom Fukui and dual descriptor indices.
 *
 * @note The quality of the resulting Fukui and dual descriptor indices heavily depends on the
 * quality of the atomic charges. We recommend Hirshfeld charges.
 *
 * @param atomicCharges The atomic charges of the structure of interest \f$q(R, N)\f$.
 * @param atomicChargesPlus The atomic charges for the same geometry with one additional electron \f$q(R, N+1)\f$.
 * @param atomicChargesMinus The atomic charges for the same geometry with one electron less \f$q(R, N-1)\f$.
 * @return A LocalConceptualDftContainer with local cDFT properties.
 */
LocalConceptualDftContainer calculateLocal(const Eigen::VectorXd& atomicCharges, const Eigen::VectorXd& atomicChargesPlus,
                                           const Eigen::VectorXd& atomicChargesMinus);

/**
 * @brief Calculates the conceptual DFT chemical potential.
 *
 * \f$\mu=\frac{E(R, N+1) - E(R, N-1)}{2}\f$
 *
 * Parr, R. G.; Pearson, R. G., J. Am. Chem. Soc. 1983, 105 (26), 7512–7516. https://doi.org/10.1021/ja00364a005.
 *
 * @param energy The energy of the structure of interest \f$E(R, N)\f$.
 * @param energyPlus The energy for the same geometry with one additional electron \f$E(R, N+1)\f$
 * @param energyMinus The energy for the same geometry with one electron less \f$E(R, N-1)\f$
 * @return The chemical potential.
 */
double calculateChemicalPotential(const double energy, const double energyPlus, const double energyMinus);

/**
 * @brief Calculates the Mulliken electronegativity.
 *
 * \f$\chi=\frac{E(R, N-1) - E(R, N+1)}{2}\f$
 *
 * Mulliken, R. S.; J. Chem. Phys. 1934, 2 (11), 782–793. https://doi.org/10.1063/1.1749394.
 *
 * @param energy The energy of the structure of interest \f$E(R, N)\f$.
 * @param energyPlus The energy for the same geometry with one additional electron \f$E(R, N+1)\f$
 * @param energyMinus The energy for the same geometry with one electron less \f$E(R, N-1)\f$
 * @return The Mulliken electronegativity.
 */
double calculateElectronegativity(const double energy, const double energyPlus, const double energyMinus);

/**
 * @brief Calculates the conceptual DFT chemical hardness.
 *
 * \f$\eta=E(R, N-1)+E(R, N+1)-2E(R, N)\f$
 *
 * Parr, R. G.; Pearson, R. G.; J. Am. Chem. Soc. 1983, 105 (26), 7512–7516. https://doi.org/10.1021/ja00364a005.
 *
 * @note In the original paper there was an extra factor of \f$\frac{1}{2}\f$ that is usually omitted today.
 * @param energy The energy of the structure of interest \f$E(R, N)\f$.
 * @param energyPlus The energy for the same geometry with one additional electron \f$E(R, N+1)\f$
 * @param energyMinus The energy for the same geometry with one electron less \f$E(R, N-1)\f$
 * @return The hardness.
 */
double calculateHardness(const double energy, const double energyPlus, const double energyMinus);

/**
 * @brief Calculates the conceptual DFT softness.
 *
 * \f$S=\frac{1}{\eta}\f$
 *
 * Yang, W.; Parr, R. G.; PNAS 1985, 82 (20), 6723–6726. https://doi.org/10.1073/pnas.82.20.6723.
 *
 * @param energy The energy of the structure of interest \f$E(R, N)\f$.
 * @param energyPlus The energy for the same geometry with one additional electron \f$E(R, N+1)\f$
 * @param energyMinus The energy for the same geometry with one electron less \f$E(R, N-1)\f$
 * @return The softness.
 */
double calculateSoftness(const double energy, const double energyPlus, const double energyMinus);

/**
 * @brief Calculates the conceptual DFT electrophilicity.
 *
 * \f$\omega=\frac{\mu^2}{2\eta}\f$
 *
 * Parr, R. G.; Szentpály, L. v.; Liu, S.; J. Am. Chem. Soc. 1999, 121 (9), 1922–1924. https://doi.org/10.1021/ja983494x.
 *
 * @param energy The energy of the structure of interest \f$E(R, N)\f$.
 * @param energyPlus The energy for the same geometry with one additional electron \f$E(R, N+1)\f$
 * @param energyMinus The energy for the same geometry with one electron less \f$E(R, N-1)\f$
 * @return The electrophilicity.
 */
double calculateElectrophilicity(const double energy, const double energyPlus, const double energyMinus);

/**
 * @brief Calculates the Fukui indices for nucleophilic attack susceptibility \f$f^+\f$.
 *
 * \f$f^+_k=q_k(R, N) - q_k(R, N+1)\f$
 *
 * Parr, R. G.; Yang, W.; J. Am. Chem. Soc. 1984, 106 (14), 4049–4050. https://doi.org/10.1021/ja00326a036.
 *
 * @note The quality of the results heavily depends on the quality of the atomic charges. We recommend Hirshfeld
 * charges.
 *
 * @param atomicCharges The atomic charges of the structure of interest \f$q(R, N)\f$.
 * @param atomicChargesPlus The atomic charges for the same geometry with one additional electron \f$q(R, N+1)\f$.
 * @param atomicChargesMinus The atomic charges for the same geometry with one electron less \f$q(R, N-1)\f$.
 * @return The Fukui indices for nucleophilic attack susceptibility.
 */
Eigen::VectorXd calculateFukuiPlus(const Eigen::VectorXd& atomicCharges, const Eigen::VectorXd& atomicChargesPlus,
                                   const Eigen::VectorXd& atomicChargesMinus);

/**
 * @brief Calculates the Fukui indices for electrophilic attack susceptibility \f$f^-\f$.
 *
 * \f$f^-_k=q_k(R, N-1) - q_k(R, N)\f$
 *
 * Parr, R. G.; Yang, W.; J. Am. Chem. Soc. 1984, 106 (14), 4049–4050. https://doi.org/10.1021/ja00326a036.
 *
 * @note The quality of the results heavily depends on the quality of the atomic charges. We recommend Hirshfeld
 * charges.
 *
 * @param atomicCharges The atomic charges of the structure of interest \f$q(R, N)\f$.
 * @param atomicChargesPlus The atomic charges for the same geometry with one additional electron \f$q(R, N+1)\f$.
 * @param atomicChargesMinus The atomic charges for the same geometry with one electron less \f$q(R, N-1)\f$.
 * @return The Fukui indices for electrophilic attack susceptibility.
 */
Eigen::VectorXd calculateFukuiMinus(const Eigen::VectorXd& atomicCharges, const Eigen::VectorXd& atomicChargesPlus,
                                    const Eigen::VectorXd& atomicChargesMinus);

/**
 * @brief Calculates the condensed to atom radical Fukui indices \f$f^0\f$.
 *
 * \f$f^0_k=\frac{1}{2}(q_k(R, N-1) - q_k(R, N+1))\f$
 *
 * Parr, R. G.; Yang, W.; J. Am. Chem. Soc. 1984, 106 (14), 4049–4050. https://doi.org/10.1021/ja00326a036.
 *
 * @note The quality of the results heavily depends on the quality of the atomic charges. We recommend Hirshfeld
 * charges.
 * @note The relevance of the radical Fukui function is known to be limited. It is included here for completeness
 * only.
 *
 * @param atomicCharges The atomic charges of the structure of interest \f$q(R, N)\f$.
 * @param atomicChargesPlus The atomic charges for the same geometry with one additional electron \f$q(R, N+1)\f$.
 * @param atomicChargesMinus The atomic charges for the same geometry with one electron less \f$q(R, N-1)\f$.
 * @return The Fukui indices for radical attacks.
 */
Eigen::VectorXd calculateFukuiRadical(const Eigen::VectorXd& atomicCharges, const Eigen::VectorXd& atomicChargesPlus,
                                      const Eigen::VectorXd& atomicChargesMinus);

/**
 * @brief Calculates the condensed to atom dual descriptor.
 *
 * \f$\Delta f_k = f_k^+ - f_k^-\f$
 *
 * Morell, C.; Grand, A.; Toro-Labbé, A.; J. Phys. Chem. A 2005, 109 (1), 205–212. https://doi.org/10.1021/jp046577a.
 *
 * @note The quality of the results heavily depends on the quality of the atomic charges. We recommend Hirshfeld
 * charges.
 *
 * @param atomicCharges The atomic charges of the structure of interest \f$q(R, N)\f$.
 * @param atomicChargesPlus The atomic charges for the same geometry with one additional electron \f$q(R, N+1)\f$.
 * @param atomicChargesMinus The atomic charges for the same geometry with one electron less \f$q(R, N-1)\f$.
 * @return The condensed to atom dual descriptor.
 */
Eigen::VectorXd calculateDualDescriptor(const Eigen::VectorXd& atomicCharges, const Eigen::VectorXd& atomicChargesPlus,
                                        const Eigen::VectorXd& atomicChargesMinus);
} // namespace ConceptualDft
} // namespace Utils
} // namespace Scine

#endif // UTILS_CONCEPTUALDFT_H
