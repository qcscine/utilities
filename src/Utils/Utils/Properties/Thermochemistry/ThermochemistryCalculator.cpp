/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ThermochemistryCalculator.h"
#include <Utils/GeometricDerivatives/NormalModeAnalysis.h>

namespace Scine {
namespace Utils {

namespace {
constexpr const double R = Constants::molarGasConstant * 1e-3 * Constants::hartree_per_kJPerMol; // hartree/K
// h * c / k_B in SI units (unit: cm * K, CODATA14)
constexpr const double wavenumberConversionCoefficient =
    Constants::planckConstant * Constants::speedOfLight / Constants::boltzmannConstant * 100;
// ln(pi * (8*pi^2*c/h)^3), the unit inside the log is m^-3 kg^-3.
// Rest of the ln is in the calculation of the rotational entropy,
// with unit length^3 mass^3, so that with the sum property of logarithms
// log(a) + log(b) = log(a*b), a*b is dimensionless.
// In this case, a is this constant, b is
//
// I_x * I_y * I_z / (\sigma^2 * (wavenumberConversionCoefficient/temperature)^3)
// |mass^3 length^6|  |----------------------- length^-3 -----------------------|
// |-------------------------- mass^3 length^3 ---------------------------------|
constexpr const double constantPartForRotational = 23.498533603003565418;
} // namespace

ThermochemistryCalculator::ThermochemistryCalculator(NormalModesContainer normalModesContainer,
                                                     Geometry::Properties::PrincipalMomentsOfInertia principalMomentsOfInertia,
                                                     ElementTypeCollection elements, int spinMultiplicity,
                                                     double electronicEnergy)
  : principalMomentsOfInertia_(std::move(principalMomentsOfInertia)),
    elements_(std::move(elements)),
    spinMultiplicity_(spinMultiplicity),
    electronicEnergy_(electronicEnergy),
    normalModesContainer_(std::move(normalModesContainer)) {
}

ThermochemistryCalculator::ThermochemistryCalculator(const HessianMatrix& hessian, const AtomCollection& atoms,
                                                     int spinMultiplicity, double electronicEnergy)
  : ThermochemistryCalculator(hessian, atoms.getElements(), atoms.getPositions(), spinMultiplicity, electronicEnergy) {
}

ThermochemistryCalculator::ThermochemistryCalculator(const HessianMatrix& hessian, ElementTypeCollection elements,
                                                     const PositionCollection& positions, int spinMultiplicity,
                                                     double electronicEnergy)
  : elements_(std::move(elements)), spinMultiplicity_(spinMultiplicity), electronicEnergy_(electronicEnergy) {
  auto masses = Utils::Geometry::Properties::getMasses(elements_);
  auto centerOfMass = Utils::Geometry::Properties::getCenterOfMass(positions, masses);
  principalMomentsOfInertia_ = Utils::Geometry::Properties::calculatePrincipalMoments(positions, masses, centerOfMass);
  normalModesContainer_ = Utils::NormalModeAnalysis::calculateNormalModes(hessian, elements_, positions);
}

ThermochemistryCalculator::ThermochemistryCalculator(ElementTypeCollection elements) : elements_(std::move(elements)) {
}

void ThermochemistryCalculator::setTemperature(double temperature) {
  if (temperature < 0.0) {
    throw std::runtime_error("A negative temperature was detected.");
  }
  temperature_ = temperature;
}

void ThermochemistryCalculator::setPressure(double pressure) {
  if (pressure < 1e-6) {
    throw std::runtime_error("Only pressures of at least 1e-6 Pa are supported.");
  }
  pressure_ = pressure;
}

void ThermochemistryCalculator::setZPVEInclusion(ZPVEInclusion inclusion) {
  zpveIncluded = inclusion;
}

ThermochemicalComponentsContainer ThermochemistryCalculator::calculate() {
  // TODO Automatically determine symmetry for all molecules
  // This only affects the symmetry of diatomic molecules
  calculateSigmaForDiatomicMolecule();
  ThermochemicalComponentsContainer container;
  container.vibrationalComponent = calculateVibrationalPart(temperature_);
  container.rotationalComponent = calculateRotationalPart(temperature_);
  container.translationalComponent = calculateTranslationalPart(temperature_, pressure_);
  container.electronicComponent = calculateElectronicPart(temperature_);
  container.overall = container.vibrationalComponent + container.rotationalComponent +
                      container.translationalComponent + container.electronicComponent;
  container.overall.symmetryNumber = sigma_;
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
  const auto& wavenumbers = getWavenumbers();

  for (auto wn : wavenumbers) {
    if (wn > 0) {
      double characteristicTemperature = wavenumberConversionCoefficient * wn;
      vibrationalTC.zeroPointVibrationalEnergy += characteristicTemperature;
      if (temperature > 1e-6) {
        double temperatureRatio = characteristicTemperature / temperature;
        double E_w = std::exp(-temperatureRatio);
        double stateProbability = 1 / (std::exp(temperatureRatio) - 1);
        vibrationalTC.enthalpy += characteristicTemperature * stateProbability;
        vibrationalTC.entropy += temperatureRatio * stateProbability - std::log(1 - E_w);
        vibrationalTC.heatCapacityP += E_w * (temperatureRatio / (E_w - 1)) * (temperatureRatio / (E_w - 1));
      }
      else {
        // At T = 0, the heat capacity vanishes, the enthalpy contributions should be only the ZPE, and the entropy
        // should vanish.
        vibrationalTC.heatCapacityP = 0;
        vibrationalTC.enthalpy = 0;
        vibrationalTC.entropy = 0;
      }
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
  vibrationalTC.gibbsFreeEnergy = vibrationalTC.enthalpy - temperature * vibrationalTC.entropy;
  return vibrationalTC;
}

ThermochemicalContainer ThermochemistryCalculator::calculateRotationalPart(double temperature) const {
  ThermochemicalContainer rotationalTC{};

  const double thermalCoefficient = wavenumberConversionCoefficient * 1e-2 / temperature * Constants::bohr_per_meter; // unit: bohr

  bool isLinear = getWavenumbers().size() == 3 * elements_.size() - 5;

  if (elements_.size() == 1 || elements_.empty()) {
    return rotationalTC;
  }
  if (isLinear) {
    double rotationalConstant = principalMomentsOfInertia_.eigenvalues(2) * Constants::electronRestMass_per_u;
    rotationalTC.enthalpy = R * temperature;
    rotationalTC.heatCapacityP = R;
    rotationalTC.heatCapacityV = rotationalTC.heatCapacityP * 3. / 5.;
    rotationalTC.entropy = R * (std::log(4.0 * Constants::pi * rotationalConstant *
                                         Constants::inverseFineStructureConstant / (sigma_ * thermalCoefficient)) +
                                1.0);
  }
  else {
    // Convert from amu * bohr^2 to m_e * bohr^2
    Eigen::Vector3d rotConst = principalMomentsOfInertia_.eigenvalues.array() * Constants::electronRestMass_per_u;

    rotationalTC.enthalpy = 1.5 * R * temperature;
    rotationalTC.heatCapacityP = 1.5 * R;
    rotationalTC.heatCapacityV = rotationalTC.heatCapacityP * 3. / 5.;
    rotationalTC.entropy =
        0.5 * R *
        (constantPartForRotational +
         std::log(rotConst(0) * rotConst(1) * rotConst(2) / (sigma_ * sigma_ * std::pow(thermalCoefficient, 3))) + 3.0);
  }
  rotationalTC.gibbsFreeEnergy = rotationalTC.enthalpy - temperature * rotationalTC.entropy;

  return rotationalTC;
}

ThermochemicalContainer ThermochemistryCalculator::calculateTranslationalPart(double temperature, double pressure) const {
  ThermochemicalContainer translationalTC{};
  double molecularMass = 0.;
  for (auto mass : Geometry::Properties::getMasses(elements_)) {
    molecularMass += mass;
  }
  translationalTC.enthalpy = 2.5 * temperature * R;
  const double pressureAtomicUnits = pressure * Constants::hartree_per_joule /
                                     (Constants::bohr_per_meter * Constants::bohr_per_meter * Constants::bohr_per_meter);
  const double massPerParticle = molecularMass * 1e-3 / Constants::avogadroNumber * 1.0 / Constants::electronRestMass;
  const double log2Pi_3over2 = std::log(2.0 * M_PI) * 3.0 / 2.0;
  const double logKb_5over2 = std::log(R) * 5.0 / 2.0;
  const double logT_5over2 = std::log(temperature) * 5.0 / 2.0;
  // q = V * (m k T /(2pi))^3/2 ; V = k T / p (atomic units: h/(2pi) = 1)
  const double logQ =
      logKb_5over2 + logT_5over2 - log2Pi_3over2 - std::log(pressureAtomicUnits) + 3.0 / 2.0 * std::log(massPerParticle);
  translationalTC.entropy = R * (5.0 / 2.0 + logQ);
  translationalTC.heatCapacityP = 2.5 * R;
  translationalTC.heatCapacityV = translationalTC.heatCapacityP * 3. / 5.;
  translationalTC.gibbsFreeEnergy = translationalTC.enthalpy - temperature * translationalTC.entropy;
  return translationalTC;
}

ThermochemicalContainer ThermochemistryCalculator::calculateElectronicPart(double temperature) const {
  ThermochemicalContainer electronicTC{};
  electronicTC.enthalpy = electronicEnergy_;
  electronicTC.entropy = R * std::log(spinMultiplicity_);
  electronicTC.gibbsFreeEnergy = electronicTC.enthalpy - temperature * electronicTC.entropy;
  return electronicTC;
}

void ThermochemistryCalculator::calculateSigmaForDiatomicMolecule() {
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
