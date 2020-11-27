/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ConceptualDft.h"

namespace Scine {
namespace Utils {

namespace ConceptualDft {

ConceptualDftContainer calculate(const double energy, const Eigen::VectorXd& atomicCharges, const double energyPlus,
                                 const Eigen::VectorXd& atomicChargesPlus, const double energyMinus,
                                 const Eigen::VectorXd& atomicChargesMinus) {
  ConceptualDftContainer container;
  // Get global components
  container.global = calculateGlobal(energy, energyPlus, energyMinus);
  // Get local descriptors
  container.local = calculateLocal(atomicCharges, atomicChargesPlus, atomicChargesMinus);

  return container;
};

GlobalConceptualDftContainer calculateGlobal(const double energy, const double energyPlus, const double energyMinus) {
  GlobalConceptualDftContainer container;
  container.chemicalPotential = calculateChemicalPotential(energy, energyPlus, energyMinus);
  container.electronegativity = calculateElectronegativity(energy, energyPlus, energyMinus);
  container.hardness = calculateHardness(energy, energyPlus, energyMinus);
  container.softness = calculateSoftness(energy, energyPlus, energyMinus);
  container.electrophilicity = calculateElectrophilicity(energy, energyPlus, energyMinus);

  return container;
};

LocalConceptualDftContainer calculateLocal(const Eigen::VectorXd& atomicCharges, const Eigen::VectorXd& atomicChargesPlus,
                                           const Eigen::VectorXd& atomicChargesMinus) {
  LocalConceptualDftContainer container;
  container.fukuiPlus = calculateFukuiPlus(atomicCharges, atomicChargesPlus, atomicChargesMinus);
  container.fukuiMinus = calculateFukuiMinus(atomicCharges, atomicChargesPlus, atomicChargesMinus);
  container.fukuiRadical = calculateFukuiRadical(atomicCharges, atomicChargesPlus, atomicChargesMinus);
  container.dualDescriptor = calculateDualDescriptor(atomicCharges, atomicChargesPlus, atomicChargesMinus);

  return container;
}

double calculateChemicalPotential(const double energy, const double energyPlus, const double energyMinus) {
  return (energyPlus - energyMinus) / 2;
};

double calculateElectronegativity(const double energy, const double energyPlus, const double energyMinus) {
  return (energyMinus - energyPlus) / 2;
};

double calculateHardness(const double energy, const double energyPlus, const double energyMinus) {
  return (energyMinus + energyPlus - 2 * energy);
};

double calculateSoftness(const double energy, const double energyPlus, const double energyMinus) {
  return 1 / (energyMinus + energyPlus - 2 * energy);
};

double calculateElectrophilicity(const double energy, const double energyPlus, const double energyMinus) {
  double chemicalPotential = calculateChemicalPotential(energy, energyPlus, energyMinus);
  double hardness = calculateHardness(energy, energyPlus, energyMinus);
  return chemicalPotential * chemicalPotential / (2 * hardness);
};

Eigen::VectorXd calculateFukuiPlus(const Eigen::VectorXd& atomicCharges, const Eigen::VectorXd& atomicChargesPlus,
                                   const Eigen::VectorXd& atomicChargesMinus) {
  Eigen::VectorXd fukuiPlus = atomicCharges - atomicChargesPlus;
  return fukuiPlus;
};

Eigen::VectorXd calculateFukuiMinus(const Eigen::VectorXd& atomicCharges, const Eigen::VectorXd& atomicChargesPlus,
                                    const Eigen::VectorXd& atomicChargesMinus) {
  Eigen::VectorXd fukuiMinus = atomicChargesMinus - atomicCharges;
  return fukuiMinus;
};

Eigen::VectorXd calculateFukuiRadical(const Eigen::VectorXd& atomicCharges, const Eigen::VectorXd& atomicChargesPlus,
                                      const Eigen::VectorXd& atomicChargesMinus) {
  Eigen::VectorXd fukuiRadical = (atomicChargesMinus - atomicChargesPlus) / 2;
  return fukuiRadical;
};

Eigen::VectorXd calculateDualDescriptor(const Eigen::VectorXd& atomicCharges, const Eigen::VectorXd& atomicChargesPlus,
                                        const Eigen::VectorXd& atomicChargesMinus) {
  Eigen::VectorXd dualDescriptor = 2 * atomicCharges - atomicChargesPlus - atomicChargesMinus;
  return dualDescriptor;
};

} // namespace ConceptualDft
} // namespace Utils
} // namespace Scine
