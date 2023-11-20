/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/CalculatorBasics/TestCalculator.h"
#include "Utils/Bonds/BondDetector.h"
#include "Utils/CalculatorBasics/PropertyList.h"
#include "Utils/CalculatorBasics/Results.h"
#include "Utils/Constants.h"
#include "Utils/GeometricDerivatives/NumericalHessianCalculator.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Settings.h"
#include "Utils/UniversalSettings/SettingsNames.h"

namespace Scine {
namespace Utils {

class TestSettings : public Settings {
 public:
  TestSettings() : Settings("TestSettings") {
    UniversalSettings::DoubleDescriptor convergence_threshold("Energy convergence limit.");
    convergence_threshold.setDefaultValue(1e-12);
    this->_fields.push_back(SettingsNames::selfConsistenceCriterion, convergence_threshold);

    UniversalSettings::IntDescriptor multiplicity("multiplicity");
    multiplicity.setDefaultValue(1);
    this->_fields.push_back(SettingsNames::spinMultiplicity, multiplicity);

    UniversalSettings::StringDescriptor spin_mode("spin mode");
    spin_mode.setDefaultValue("restricted");
    this->_fields.push_back(SettingsNames::spinMode, spin_mode);

    this->resetToDefaults();
  };
  ~TestSettings() override = default;
};

TestCalculator::TestCalculator() {
  _settings = std::make_shared<TestSettings>();
}
bool TestCalculator::supportsMethodFamily(const std::string& methodFamily) const {
  return methodFamily == "TEST";
}
void TestCalculator::setStructure(const AtomCollection& structure) {
  _structure = structure;
  _results = Results{};
}
void TestCalculator::modifyPositions(PositionCollection newPositions) {
  _structure.setPositions(newPositions);
  _results = Results{};
}
const PositionCollection& TestCalculator::getPositions() const {
  return _structure.getPositions();
}
void TestCalculator::setRequiredProperties(const PropertyList& requiredProperties) {
  _requiredProperties = requiredProperties;
}
PropertyList TestCalculator::getRequiredProperties() const {
  return _requiredProperties;
}
PropertyList TestCalculator::possibleProperties() const {
  return Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::Hessian | Utils::Property::BondOrderMatrix;
}
const Results& TestCalculator::calculate(std::string /*dummy*/) {
  auto positions = _structure.getPositions();
  auto elements = _structure.getElements();
  const unsigned int nAtoms = _structure.size();
  double energy = 0.0;
  GradientCollection g(_structure.size(), 3);
  g.setZero();

  for (unsigned int i = 0; i < nAtoms; i++) {
    auto posI = positions.row(i);
    const double radi = ElementInfo::covalentRadius(elements[i]);
    for (unsigned int j = 0; j < i; j++) {
      /*
       * Energy
       */
      const Eigen::Vector3d r = posI - positions.row(j);
      const double dist = r.norm();
      const double radj = ElementInfo::covalentRadius(elements[j]);
      const double rMin = radi + radj;
      const double lj = rMin / dist;
      const double lj6 = lj * lj * lj * lj * lj * lj;
      const double lj12 = lj6 * lj6;
      // Scaling factor for all pairs, a common c-c bond is roughly 1.0
      // max scaling 2 in order to not have unbreakable bonds for heavy atoms
      const double scaling = std::min(rMin / 2.0, 2.0);
      // Gaussian exponent
      const double gauss = (dist - 2.5 * scaling) / (1.0 * scaling);
      // Gaussian exponent evaluated and scaled
      const double egauss = (0.4 / dist) * scaling * std::exp(-gauss * gauss);
      const double wellDepth = 0.2 * scaling;
      // compose energy contribution
      energy += wellDepth * (lj12 - 2.0 * lj6) + egauss;
      // truncate energy
      energy = truncateOff(energy);
      /*
       * Gradient
       */
      // Gaussian derivative (all but the Gaussian it self)
      const double gderiv = -1.0 * (-5.0 * scaling * dist + 2.0 * dist * dist + 1) / dist;
      // Partial derivative of E w.r.t. r
      const double rDeriv = gderiv * egauss + wellDepth * 12.0 * (lj6 / dist - lj12 / dist);
      // Derivative of E w.r.t. x,y,z
      g(i, 0) += (rDeriv / dist) * r[0];
      g(i, 1) += (rDeriv / dist) * r[1];
      g(i, 2) += (rDeriv / dist) * r[2];
      g(j, 0) -= (rDeriv / dist) * r[0];
      g(j, 1) -= (rDeriv / dist) * r[1];
      g(j, 2) -= (rDeriv / dist) * r[2];
    }
  }
  // truncate gradients
  for (unsigned int i = 0; i < nAtoms; i++) {
    g(i, 0) = truncateOff(g(i, 0));
    g(i, 1) = truncateOff(g(i, 1));
    g(i, 2) = truncateOff(g(i, 2));
  }
  _results = Results();
  _results.set<Property::SuccessfulCalculation>(true);
  _results.set<Property::Energy>(energy);
  if (_settings->getInt(SettingsNames::spinMultiplicity) != 1) {
    _results.set<Property::Energy>(energy - _settings->getInt(SettingsNames::spinMultiplicity));
  }
  _results.set<Property::Gradients>(g);
  if (_requiredProperties.containsSubSet(Scine::Utils::Property::BondOrderMatrix)) {
    auto bos = BondDetector::detectBonds(_structure);
    _results.set<Property::BondOrderMatrix>(bos);
  }
  if (_requiredProperties.containsSubSet(Scine::Utils::Property::Hessian)) {
    auto copy(*this);
    NumericalHessianCalculator nhCalc(copy);
    auto hessResults = nhCalc.calculate();
    _results.set<Property::Hessian>(hessResults.get<Property::Hessian>());
  }
  return _results;
}
std::string TestCalculator::name() const {
  return std::string("TestCalculator");
}
const Settings& TestCalculator::settings() const {
  return *_settings;
}
Settings& TestCalculator::settings() {
  return *_settings;
}
std::shared_ptr<Core::State> TestCalculator::getState() const {
  return nullptr;
}
void TestCalculator::loadState(std::shared_ptr<Core::State> /* state */) {
}
Utils::Results& TestCalculator::results() {
  return _results;
}
const Utils::Results& TestCalculator::results() const {
  return _results;
}
std::unique_ptr<Utils::AtomCollection> TestCalculator::getStructure() const {
  return std::make_unique<AtomCollection>(_structure);
}
void TestCalculator::generateWavefunctionInformation(const std::string& out) {
  std::ofstream fileOut(out);
  if (!fileOut.is_open()) {
    throw std::runtime_error("Impossible to open file " + out + ".");
  }
  generateWavefunctionInformation(fileOut);
}
void TestCalculator::generateWavefunctionInformation(std::ostream& out) {
  out << "This is a test wavefunction output information." << std::endl;
}
void TestCalculator::setPrecision(double value) {
  _precision = value;
}
double TestCalculator::getPrecision() {
  return _precision;
}
double TestCalculator::truncateOff(double value) {
  double factor_10 = std::pow(10.0, _precision);
  return std::trunc(value * factor_10) / factor_10;
}

} // namespace Utils
} // namespace Scine
