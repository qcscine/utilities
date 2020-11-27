/**
 * @file ThermochemistryContainer.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_THERMOCHEMISTRYCALCULATOR_H
#define UTILS_THERMOCHEMISTRYCALCULATOR_H

#include <Utils/GeometricDerivatives/NormalModesContainer.h>
#include <Utils/Geometry.h>
namespace Scine {
namespace Utils {

/**
 * @brief Struct containing the single thermochemical properties of interest.
 */
struct ThermochemicalContainer {
  double entropy = 0, enthalpy = 0, heatCapacityP = 0, heatCapacityV = 0, gibbsFreeEnergy = 0, zeroPointVibrationalEnergy = 0;
  int symmetryNumber = 1;
  ThermochemicalContainer operator+(const ThermochemicalContainer& rhs) const {
    ThermochemicalContainer result = *this;
    result.symmetryNumber = this->symmetryNumber;
    result.entropy += rhs.entropy;
    result.enthalpy += rhs.enthalpy;
    result.heatCapacityP += rhs.heatCapacityP;
    result.heatCapacityV += rhs.heatCapacityV;
    result.gibbsFreeEnergy += rhs.gibbsFreeEnergy;
    result.zeroPointVibrationalEnergy += rhs.zeroPointVibrationalEnergy;
    return result;
  }
  /** @brief Overloaded multiplication for unit conversion. */
  ThermochemicalContainer operator*(double factor) const {
    ThermochemicalContainer result = *this;
    result.entropy *= factor;
    result.enthalpy *= factor;
    result.heatCapacityP *= factor;
    result.heatCapacityV *= factor;
    result.gibbsFreeEnergy *= factor;
    result.zeroPointVibrationalEnergy *= factor;
    return result;
  }
};

/**
 * @brief Struct containing the vibrational, rotational, translational and the overall thermochemical components.
 */
struct ThermochemicalComponentsContainer {
  ThermochemicalContainer vibrationalComponent{}, rotationalComponent{}, translationalComponent{},
      electronicComponent{}, overall{};
  /** @brief Overloaded multiplication for unit conversion */
  ThermochemicalComponentsContainer operator*(double factor) const {
    ThermochemicalComponentsContainer result;
    result.vibrationalComponent = this->vibrationalComponent * factor;
    result.rotationalComponent = this->rotationalComponent * factor;
    result.translationalComponent = this->translationalComponent * factor;
    result.electronicComponent = this->electronicComponent * factor;
    result.overall = this->overall * factor;
    return result;
  }
};

/**
 * In NDDO semiempirical methods, the ZPVE is already included in the electronic energy
 * due to the way they are parametrized. In the thermochemical calculation it must then
 * not be included again.
 */
enum class ZPVEInclusion { alreadyIncluded, notIncluded };

/**
 * @class ThermochemistryCalculator @file ThermochemistryCalculator.h
 * @brief This class calculates and stores the most important thermochemical data.
 * The class calculates important thermochemical descriptors, as the Zero Point Vibrational Energy,
 * the standard enthalpy of formation, the system entropy, the heat capacity and the Gibbs' free energy
 * from the results of a vibrational analysis.
 *
 * By default a temperature of 298.15 K and a molecular symmetry number of one are assumed. They can be adapted via
 * the functions setTemperature() and setMolecularSymmetryNumber(), respectively. Using an incorrect symmetry number
 * results into a wrong rotational entropy.
 */
class ThermochemistryCalculator {
 public:
  explicit ThermochemistryCalculator(ElementTypeCollection elements);
  ThermochemistryCalculator(NormalModesContainer normalModesContainer,
                            Geometry::PrincipalMomentsOfInertia principalMomentsOfInertia,
                            ElementTypeCollection elements, int spinMulitplicity, double electronicEnergy);
  ThermochemistryCalculator(const HessianMatrix& hessian, ElementTypeCollection elements,
                            const PositionCollection& positions, int spinMultiplicity, double electronicEnergy);
  ~ThermochemistryCalculator() = default;

  /**
   * @brief Setter for the temperature at which the thermochemical calculation is performed.
   */
  void setTemperature(double temperature);
  /**
   * In NDDO semiempirical methods, the ZPVE is already included in the electronic energy
   * due to the way they are parametrized. In the thermochemical calculation it must then
   * not be included again.
   * @param inclusion notIncluded is the standard way, alreadyIncluded is in case of the NDDO semiempirical methods.
   */
  void setZPVEInclusion(ZPVEInclusion inclusion);

  ThermochemicalComponentsContainer calculate();

  /**
   * @brief Sets the symmetry sigma factor.
   * @param sigma The symmetry factor related to the point group symmetry of the molecule.
   * Examples:
   * Cn,v/h : n
   * Dn,v/h : 2*n
   * C_inf,v : 1
   * C_inf,h : 2
   * S_n : n/2
   * T : 12
   * O : 24
   * I : 60
   */
  void setMolecularSymmetryNumber(int sigma);

 protected:
  std::vector<double> getWavenumbers() const;
  Geometry::PrincipalMomentsOfInertia principalMomentsOfInertia_;
  ElementTypeCollection elements_;
  double temperature_{298.15};
  int spinMultiplicity_{1};
  double electronicEnergy_{0};
  int sigma_{1};
  ZPVEInclusion zpveIncluded{ZPVEInclusion::notIncluded};

 private:
  ThermochemicalContainer calculateVibrationalPart(double temperature) const;
  ThermochemicalContainer calculateRotationalPart(double temperature) const;
  ThermochemicalContainer calculateTranslationalPart(double temperature) const;
  ThermochemicalContainer calculateElectronicPart(double temperature) const;
  // Sets the symmetry number for diatomic molecules
  void calculateSigmaForDiatomicMolecule();
  // TODO Automatically determine symmetry for all molecules
  NormalModesContainer normalModesContainer_{};
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_THERMOCHEMISTRYCALCULATOR_H
