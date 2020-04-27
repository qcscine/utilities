/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_LCAOMETHOD_H
#define UTILS_LCAOMETHOD_H

#include "SinglePointMethod.h"
#include <Utils/DataStructures/AtomsOrbitalsIndexes.h>
#include <Utils/DataStructures/DensityMatrix.h>
#include <Utils/DataStructures/MolecularOrbitals.h>
#include <Utils/DataStructures/SingleParticleEnergies.h>
#include <Utils/DataStructures/SpinAdaptedMatrix.h>
#include <Utils/Math/DerivOrderEnum.h>
#include <Utils/Scf/LcaoUtils/ElectronicOccupation.h>
#include <Eigen/Core>
#include <memory>

namespace Scine {
namespace Utils {

namespace LcaoUtils {
class ElectronicOccupationGenerator;
}
class ElectronicContributionCalculator;
class OverlapCalculator;
class RepulsionCalculator;
class StructureDependentInitializer;
class AdditiveElectronicContribution;

class UnrestrictedCalculationNotAvailableException : public std::exception {
 public:
  const char* what() const noexcept final {
    return "Unrestricted calculation not possible with this method.";
  }
};

/*!
 * Base class describing a method based on a LCAO ansatz.
 */
class LcaoMethod : public SinglePointMethod {
 public:
  explicit LcaoMethod(bool unrestrictedCalculationPossible, Utils::derivOrder maximalOrder, bool orthogonalBasisSet = false);
  virtual ~LcaoMethod() override;

  /*! Initialize the method <b>after</b> the parameters have been set or loaded. */
  virtual void initialize();

  void calculate(Utils::derivativeType d) override;

  void calculateDensityIndependentQuantities(Utils::derivativeType d = Utils::derivativeType::first);
  void assembleFockMatrix();
  void computeEnergyAndDerivatives(Utils::derivativeType d = Utils::derivativeType::first);

  /*! \return true if the basis set is orthogonal and the electronic structure is optimized in
      a normal eigenvalue problem, false if it is a generalized eigenvalue problem. */
  bool basisSetIsOrthogonal() const;

  /**
   * @brief True if only the eigenpairs corresponding to the occupied MOs need to be calculated.
   * Subspace calculation is then performed with the Davidson-Liu Algorithm.
   */
  bool solvesOnlyOccupiedManifold() const;
  void setOnlyOccupiedManifoldToSolve(bool onlyOccupiedToSolve);

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
  const LcaoUtils::ElectronicOccupation& getElectronicOccupation() const;

  /*! Set how the electronic occupation is to be set. By default, the Aufbau principle is used. */
  void setElectronicOccupationGenerator(std::unique_ptr<LcaoUtils::ElectronicOccupationGenerator>&& electronicOccupationSetter);

  /** @brief Adds an electronic contribution to the Fock matrix. */
  void addElectronicContribution(std::shared_ptr<AdditiveElectronicContribution> contribution);

  void verifyPesValidity();

  void calculateDensity();

 protected:
  /*! Resets the size of the different matrices. To be called during initialization of the specific method. */
  void resizeLcaoMethodMatrices();
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
  LcaoUtils::ElectronicOccupation occupation_;

  int molecularCharge_ = 0;
  bool unrestrictedCalculationPossible_;
  bool unrestrictedCalculationRunning_ = false;
  bool solveOnlyOccupiedManifold_ = false;
  int spinMultiplicity_ = 1;
  std::vector<double> coreCharges_; // Atomic core charges
  std::unique_ptr<LcaoUtils::ElectronicOccupationGenerator> electronicOccupationGenerator_;

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

inline bool LcaoMethod::basisSetIsOrthogonal() const {
  return basisSetIsOrthogonal_;
}

inline const SpinAdaptedMatrix& LcaoMethod::getFockMatrix() const {
  return fockMatrix_;
}

inline const DensityMatrix& LcaoMethod::getDensityMatrix() const {
  return densityMatrix_;
}

inline const Eigen::MatrixXd& LcaoMethod::getEnergyWeightedDensityMatrix() const {
  return energyWeightedDensityMatrix_;
}

inline MolecularOrbitals& LcaoMethod::getMolecularOrbitals() {
  return eigenvectorMatrix_;
}

inline const MolecularOrbitals& LcaoMethod::getMolecularOrbitals() const {
  return eigenvectorMatrix_;
}

inline int LcaoMethod::getNumberElectrons() const {
  return nElectrons_;
}

inline int LcaoMethod::getNumberAtomicOrbitals() const {
  return nAOs_;
}

inline void LcaoMethod::setMolecularCharge(int c) {
  molecularCharge_ = c;
}

inline int LcaoMethod::getMolecularCharge() const {
  return molecularCharge_;
}

inline int LcaoMethod::spinMultiplicity() const {
  return spinMultiplicity_;
}

inline bool LcaoMethod::unrestrictedCalculationRunning() const {
  return unrestrictedCalculationRunning_;
}

inline bool LcaoMethod::unrestrictedCalculationPossible() const {
  return unrestrictedCalculationPossible_;
}

inline const AtomsOrbitalsIndexes& LcaoMethod::getAtomsOrbitalsIndexesHolder() const {
  return aoIndexes_;
}

inline const SingleParticleEnergies& LcaoMethod::getSingleParticleEnergies() const {
  return singleParticleEnergies_;
}

inline const LcaoUtils::ElectronicOccupation& LcaoMethod::getElectronicOccupation() const {
  return occupation_;
}

inline double LcaoMethod::getRepulsionEnergy() const {
  return repulsionEnergy_;
}

inline double LcaoMethod::getElectronicEnergy() const {
  return electronicEnergy_;
}

} // namespace Utils
} // namespace Scine
#endif // UTILS_LCAOMETHOD_H
