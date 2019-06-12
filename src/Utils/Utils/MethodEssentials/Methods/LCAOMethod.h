/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_LCAOMETHOD_H
#define UTILS_LCAOMETHOD_H

#include "SinglePointMethod.h"
#include <Utils/Math/DerivOrderEnum.h>
#include <Utils/MethodEssentials/util/AtomsOrbitalsIndexes.h>
#include <Utils/MethodEssentials/util/DensityMatrix.h>
#include <Utils/MethodEssentials/util/LcaoUtil/ElectronicOccupation.h>
#include <Utils/MethodEssentials/util/MolecularOrbitals.h>
#include <Utils/MethodEssentials/util/SingleParticleEnergies.h>
#include <Utils/MethodEssentials/util/SpinAdaptedMatrix.h>
#include <Eigen/Core>
#include <memory>

namespace Scine {
namespace Utils {

namespace LcaoUtil {
class ElectronicOccupationGenerator;
}
class ElectronicContributionCalculator;
class OverlapCalculator;
class RepulsionCalculator;
class StructureDependentInitializer;

/*!
 * Base class describing a method based on a LCAO ansatz.
 */
class LCAOMethod : public SinglePointMethod {
 public:
  explicit LCAOMethod(bool unrestrictedCalculationPossible, Utils::derivOrder maximalOrder, bool orthogonalBasisSet = false);
  virtual ~LCAOMethod() override;

  /*! Initialize the method <b>after</b> the parameters have been set or loaded. */
  virtual void initialize();

  void calculate(Utils::derivativeType d) override;

  void calculateDensityIndependentQuantities(Utils::derivativeType d = Utils::derivativeType::first);
  void assembleFockMatrix();
  void computeEnergyAndDerivatives(Utils::derivativeType d = Utils::derivativeType::first);

  /*! \return true if the basis set is orthogonal and the electronic structure is optimized in
      a normal eigenvalue problem, false if it is a generalized eigenvalue problem. */
  bool basisSetIsOrthogonal() const;

  /*! \name Getters and setters for matrices
   * @{
   */
  const Eigen::MatrixXd& getOverlapMatrix() const;
  void setOverlapMatrix(const Eigen::MatrixXd& S);

  const SpinAdaptedMatrix& getFockMatrix() const;
  void setFockMatrix(SpinAdaptedMatrix F);

  const DensityMatrix& getDensityMatrix() const;
  void setDensityMatrix(DensityMatrix P);

  MolecularOrbitals& getMolecularOrbitals();
  const MolecularOrbitals& getMolecularOrbitals() const;
  void setMolecularOrbitals(MolecularOrbitals C);

  const Eigen::MatrixXd& getEnergyWeightedDensityMatrix() const;
  void setEnergyWeightedDensityMatrix(const Eigen::MatrixXd& W);
  /*! @} */

  int getNumberElectrons() const;
  int getNumberAtomicOrbitals() const;

  void setMolecularCharge(int c);
  int getMolecularCharge() const;
  void setSpinMultiplicity(int s);
  int spinMultiplicity() const;
  void setUnrestrictedCalculation(bool b);
  bool unrestrictedCalculationRunning() const;
  bool unrestrictedCalculationPossible() const;

  double getElectronicEnergy() const;
  double getRepulsionEnergy() const;

  /*! Get the HOMO-LUMO gap (Unit: Hartree). */
  double getHomoLumoGap() const;

  /*! Returns a const reference to the object containing information about which atoms have how many orbitals and what
   * their indexes are. */
  const AtomsOrbitalsIndexes& getAtomsOrbitalsIndexesHolder() const;
  /*! Returns a const reference to the object containing information about single-particle energies. */
  const SingleParticleEnergies& getSingleParticleEnergies() const;
  /*! Returns a const reference to the object containing information about the electronic occupation. */
  const LcaoUtil::ElectronicOccupation& getElectronicOccupation() const;

  /*! Set how the electronic occupation is to be set. By default, the Aufbau principle is used. */
  void setElectronicOccupationGenerator(std::unique_ptr<LcaoUtil::ElectronicOccupationGenerator>&& electronicOccupationSetter);

  void verifyPesValidity();

  void calculateDensity();

 protected:
  /*! Resets the size of the different matrices. To be called during initialization of the specific method. */
  void resizeLCAOMethodMatrices();
  void calculateOccupation();
  void calculateEnergyWeightedDensity();
  void calculateBondOrderMatrix();
  void calculateAtomicCharges();

  Eigen::MatrixXd overlapMatrix_;
  Eigen::MatrixXd energyWeightedDensityMatrix_;
  SpinAdaptedMatrix fockMatrix_;
  MolecularOrbitals eigenvectorMatrix_;
  DensityMatrix densityMatrix_;
  int nAOs_;
  int nElectrons_ = 0;
  int nElectronsForUnchargedSpecies_ = 0;
  SingleParticleEnergies singleParticleEnergies_;
  AtomsOrbitalsIndexes aoIndexes_;
  LcaoUtil::ElectronicOccupation occupation_;

  int molecularCharge_ = 0;
  bool unrestrictedCalculationPossible_;
  bool unrestrictedCalculationRunning_ = false;
  int spinMultiplicity_ = 1;
  std::vector<double> coreCharges_; // Atomic core charges
  std::unique_ptr<LcaoUtil::ElectronicOccupationGenerator> electronicOccupationGenerator_;

  std::shared_ptr<StructureDependentInitializer> initializer_;
  std::shared_ptr<RepulsionCalculator> rep_;
  std::shared_ptr<OverlapCalculator> overlapCalculator_;
  std::shared_ptr<ElectronicContributionCalculator> electronicPart_;
  double electronicEnergy_, repulsionEnergy_;

 private:
  void verifyChargeValidity();
  void verifyMultiplicityValidity();
  void verifyUnrestrictedValidity();
  void calculateOccupationAndDensity();
  template<Utils::derivativeType O>
  void calculateDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<O>& derivatives);

  const bool basisSetIsOrthogonal_;
};

inline bool LCAOMethod::basisSetIsOrthogonal() const {
  return basisSetIsOrthogonal_;
}

inline const SpinAdaptedMatrix& LCAOMethod::getFockMatrix() const {
  return fockMatrix_;
}

inline const DensityMatrix& LCAOMethod::getDensityMatrix() const {
  return densityMatrix_;
}

inline const Eigen::MatrixXd& LCAOMethod::getEnergyWeightedDensityMatrix() const {
  return energyWeightedDensityMatrix_;
}

inline MolecularOrbitals& LCAOMethod::getMolecularOrbitals() {
  return eigenvectorMatrix_;
}

inline const MolecularOrbitals& LCAOMethod::getMolecularOrbitals() const {
  return eigenvectorMatrix_;
}

inline int LCAOMethod::getNumberElectrons() const {
  return nElectrons_;
}

inline int LCAOMethod::getNumberAtomicOrbitals() const {
  return nAOs_;
}

inline void LCAOMethod::setMolecularCharge(int c) {
  molecularCharge_ = c;
}

inline int LCAOMethod::getMolecularCharge() const {
  return molecularCharge_;
}

inline int LCAOMethod::spinMultiplicity() const {
  return spinMultiplicity_;
}

inline bool LCAOMethod::unrestrictedCalculationRunning() const {
  return unrestrictedCalculationRunning_;
}

inline bool LCAOMethod::unrestrictedCalculationPossible() const {
  return unrestrictedCalculationPossible_;
}

inline const AtomsOrbitalsIndexes& LCAOMethod::getAtomsOrbitalsIndexesHolder() const {
  return aoIndexes_;
}

inline const SingleParticleEnergies& LCAOMethod::getSingleParticleEnergies() const {
  return singleParticleEnergies_;
}

inline const LcaoUtil::ElectronicOccupation& LCAOMethod::getElectronicOccupation() const {
  return occupation_;
}

inline double LCAOMethod::getRepulsionEnergy() const {
  return repulsionEnergy_;
}

inline double LCAOMethod::getElectronicEnergy() const {
  return electronicEnergy_;
}

} // namespace Utils
} // namespace Scine
#endif // UTILS_LCAOMETHOD_H
