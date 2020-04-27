/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ThermochemistryCalculator.h"
#include <Utils/GeometricDerivatives/NormalModeAnalyzer.h>

namespace Scine {
namespace Utils {

namespace {
constexpr const double R = 0.008314510 * Constants::hartree_per_kJPerMol; // hartree/K
// h * c / k_B in SI units
constexpr const double wavenumberConversionCoefficient = 6.6260755e-34 * 299792458 / 1.380658e-23 * 100;
constexpr const double invCm_per_AngstromSquaredAndAmu = 16.8576304198232;
} // namespace

ThermochemistryCalculator::ThermochemistryCalculator(NormalModesContainer normalModesContainer,
                                                     Geometry::PrincipalMomentsOfInertia principalMomentsOfInertia,
                                                     ElementTypeCollection elements, int spinMultiplicity,
                                                     double electronicEnergy)
  : normalModesContainer_(std::move(normalModesContainer)),
    principalMomentsOfInertia_(std::move(principalMomentsOfInertia)),
    elements_(std::move(elements)),
    spinMultiplicity_(spinMultiplicity),
    electronicEnergy_(electronicEnergy),
    sigma_(1) {
}

ThermochemistryCalculator::ThermochemistryCalculator(const HessianMatrix& hessian, ElementTypeCollection elements,
                                                     const PositionCollection& positions, int spinMultiplicity,
                                                     double electronicEnergy)
  : elements_(std::move(elements)), spinMultiplicity_(spinMultiplicity), electronicEnergy_(electronicEnergy), sigma_(1) {
  auto masses = Utils::Geometry::getMasses(elements_);
  auto centerOfMass = Utils::Geometry::getCenterOfMass(positions, masses);
  principalMomentsOfInertia_ = Utils::Geometry::calculatePrincipalMoments(positions, masses, centerOfMass);
  Utils::NormalModeAnalyzer analyzer(hessian, elements_, positions);
  normalModesContainer_ = analyzer.calculateNormalModes();
}

ThermochemistryCalculator::ThermochemistryCalculator(ElementTypeCollection elements) : elements_(std::move(elements)) {
}

void ThermochemistryCalculator::setTemperature(double temperature) {
  temperature_ = temperature;
}

void ThermochemistryCalculator::setZPVEInclusion(ZPVEInclusion inclusion) {
  zpveIncluded = inclusion;
}

ThermochemicalComponentsContainer ThermochemistryCalculator::calculate() {
  calculateSigmaForLinearMolecule();
  ThermochemicalComponentsContainer container;
  container.overall.temperature = temperature_;
  container.vibrationalComponent = calculateVibrationalPart(temperature_);
  container.rotationalComponent = calculateRotationalPart(temperature_);
  container.translationalComponent = calculateTranslationalPart(temperature_);
  container.electronicComponent = calculateElectronicPart(temperature_);
  container.overall = container.vibrationalComponent + container.rotationalComponent +
                      container.translationalComponent + container.electronicComponent;

  return container;
}

std::vector<double> ThermochemistryCalculator::getWavenumbers() const {
  return normalModesContainer_.getWaveNumbers();
}

void ThermochemistryCalculator::setMolecularSymmetryNumber(int sigma) {
  sigma_ = sigma;
}

ThermochemicalContainer ThermochemistryCalculator::calculateVibrationalPart(double temperature) const {
  ThermochemicalContainer vibrationalTC{};
  vibrationalTC.temperature = temperature;
  const auto& wavenumbers = getWavenumbers();

  for (auto wn : wavenumbers) {
    if (wn > 0) {
      double characteristicTemperature = wavenumberConversionCoefficient * wn;
      double temperatureRatio = characteristicTemperature / vibrationalTC.temperature;
      double E_w = std::exp(-temperatureRatio);
      double stateProbability = 1 / (std::exp(temperatureRatio) - 1);
      vibrationalTC.zeroPointVibrationalEnergy += characteristicTemperature;
      vibrationalTC.enthalpy += characteristicTemperature * stateProbability;
      vibrationalTC.entropy += temperatureRatio * stateProbability - std::log(1 - E_w);
      vibrationalTC.heatCapacityP += E_w * (temperatureRatio / (E_w - 1)) * (temperatureRatio / (E_w - 1));
    }
  }
  vibrationalTC.zeroPointVibrationalEnergy *= 0.5 * R;
  vibrationalTC.enthalpy *= R;
  if (zpveIncluded == ZPVEInclusion::notIncluded) {
    vibrationalTC.enthalpy += vibrationalTC.zeroPointVibrationalEnergy;
  }
  vibrationalTC.entropy *= R;
  vibrationalTC.heatCapacityP *= R;
  vibrationalTC.heatCapacityV = vibrationalTC.heatCapacityP * 3. / 5.;
  vibrationalTC.gibbsFreeEnergy = vibrationalTC.enthalpy - vibrationalTC.temperature * vibrationalTC.entropy;
  return vibrationalTC;
}

ThermochemicalContainer ThermochemistryCalculator::calculateRotationalPart(double temperature) const {
  ThermochemicalContainer rotationalTC{};
  rotationalTC.temperature = temperature;

  const double thermalCoefficient = wavenumberConversionCoefficient / rotationalTC.temperature;

  bool isLinear = getWavenumbers().size() == 3 * elements_.size() - 5;

  if (elements_.size() == 1 || elements_.empty()) {
    return rotationalTC;
  }
  else if (isLinear) {
    double rotationalConstant = principalMomentsOfInertia_.eigenvalues(2) * invCm_per_AngstromSquaredAndAmu;
    rotationalTC.enthalpy = R * rotationalTC.temperature;
    rotationalTC.heatCapacityP = R;
    rotationalTC.heatCapacityV = rotationalTC.heatCapacityP * 3. / 5.;
    rotationalTC.entropy = (R * std::log(1 / (sigma_ * thermalCoefficient * rotationalConstant)) + R);
  }
  else {
    Eigen::Vector3d rotConst = principalMomentsOfInertia_.eigenvalues.array() * invCm_per_AngstromSquaredAndAmu;
    rotationalTC.enthalpy = 1.5 * R * rotationalTC.temperature;
    rotationalTC.heatCapacityP = 1.5 * R;
    rotationalTC.heatCapacityV = rotationalTC.heatCapacityP * 3. / 5.;
    rotationalTC.entropy = 0.5 * R *
                           (std::log(Constants::pi / (sigma_ * sigma_ * rotConst(0) * rotConst(1) * rotConst(2) *
                                                      std::pow(thermalCoefficient, 3))) +
                            3);
  }
  rotationalTC.gibbsFreeEnergy = rotationalTC.enthalpy - rotationalTC.temperature * rotationalTC.entropy;

  return rotationalTC;
}

ThermochemicalContainer ThermochemistryCalculator::calculateTranslationalPart(double temperature) const {
  ThermochemicalContainer translationalTC{};
  double molecularMass = 0.;
  for (auto mass : Geometry::getMasses(elements_)) {
    molecularMass += mass;
  }
  // From MOPAC manual
  translationalTC.temperature = temperature;
  translationalTC.enthalpy = 2.5 * translationalTC.temperature * R;
  translationalTC.entropy =
      (9.93608e-4 * (5 * std::log(translationalTC.temperature) + 3 * std::log(molecularMass)) - 2.31482e-3) *
      Constants::hartree_per_kCalPerMol;
  translationalTC.heatCapacityP = 2.5 * R;
  translationalTC.heatCapacityV = translationalTC.heatCapacityP * 3. / 5.;
  translationalTC.gibbsFreeEnergy = translationalTC.enthalpy - translationalTC.temperature * translationalTC.entropy;
  return translationalTC;
}

ThermochemicalContainer ThermochemistryCalculator::calculateElectronicPart(double temperature) const {
  ThermochemicalContainer electronicTC{};
  electronicTC.temperature = temperature;
  electronicTC.enthalpy = electronicEnergy_;
  electronicTC.entropy = R * std::log(spinMultiplicity_);
  electronicTC.gibbsFreeEnergy = electronicTC.enthalpy - electronicTC.temperature * electronicTC.entropy;
  return electronicTC;
}

void ThermochemistryCalculator::calculateSigmaForLinearMolecule() {
  if (elements_.size() == 2) {
    if (elements_[1] == elements_[0]) {
      setMolecularSymmetryNumber(2);
    }
    else {
      setMolecularSymmetryNumber(1);
    }
  }
}

} // namespace Utils
} // namespace Scine
