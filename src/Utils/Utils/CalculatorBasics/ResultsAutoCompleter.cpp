/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/CalculatorBasics/ResultsAutoCompleter.h>
#include <Utils/Properties/Thermochemistry/ThermochemistryCalculator.h>
#include <Utils/Scf/LcaoUtils/DensityMatrixGenerator.h>
#include <Utils/Scf/LcaoUtils/ElectronicOccupation.h>
#include <Utils/Scf/LcaoUtils/LcaoUtils.h>

namespace Scine {
namespace Utils {
// Constructor
ResultsAutoCompleter::ResultsAutoCompleter(AtomCollection& atomCollection) {
  // Initialise map with lists of properties required to autocompute a certain other property
  fillRequiredPropertiesMap();
  // Set the core charges equal to the nuclear charges
  setCoreChargesToNuclearCharges(atomCollection);
  // Make the default choice for wantedProperties_
  setDefaultWantedProperties();
}

void ResultsAutoCompleter::fillRequiredPropertiesMap() {
  requiredPropertiesMap_[Property::DensityMatrix] =
      Property::SuccessfulCalculation | Property::ElectronicOccupation | Property::CoefficientMatrix;
  requiredPropertiesMap_[Property::Thermochemistry] =
      Property::SuccessfulCalculation | Property::Hessian | Property::ElectronicOccupation; // Electronic occ. for spin
  requiredPropertiesMap_[Property::AtomicCharges] =
      Property::SuccessfulCalculation | Property::DensityMatrix | Property::OverlapMatrix | Property::AOtoAtomMapping;
  requiredPropertiesMap_[Property::BondOrderMatrix] =
      Property::SuccessfulCalculation | Property::DensityMatrix | Property::OverlapMatrix | Property::AOtoAtomMapping;
}

void ResultsAutoCompleter::setDefaultWantedProperties() {
  wantedProperties_ = Property::DensityMatrix | Property::AtomicCharges | Property::BondOrderMatrix | Property::Thermochemistry;
}

void ResultsAutoCompleter::setTemperature(double temperature) {
  temperature_ = temperature;
}

void ResultsAutoCompleter::setMolecularSymmetryNumber(int sigma) {
  sigma_ = sigma;
}

void ResultsAutoCompleter::setZPVEInclusion(ZPVEInclusion zpveInclusion) {
  zpveIncluded_ = zpveInclusion;
}

void ResultsAutoCompleter::addOneWantedProperty(Scine::Utils::Property wantedProperty) {
  wantedProperties_.addProperty(wantedProperty);
}

void ResultsAutoCompleter::setWantedProperties(PropertyList wantedProperties) {
  wantedProperties_ = wantedProperties;
}

void ResultsAutoCompleter::setCoreChargesToNuclearCharges(Scine::Utils::AtomCollection& atomCollection) {
  ElementTypeCollection elementTypeCollection = atomCollection.getElements();
  coreCharges_.clear();
  for (auto const& elementType : elementTypeCollection) {
    coreCharges_.emplace_back(ElementInfo::Z(elementType));
  }
}

void ResultsAutoCompleter::setCoreCharges(std::vector<double> coreCharges) {
  coreCharges_ = std::move(coreCharges);
}

void ResultsAutoCompleter::setAllPropertiesAsWanted() {
  for (auto property : allProperties) {
    wantedProperties_.addProperty(property);
  }
}

bool ResultsAutoCompleter::propertyGeneratable(const Results& results, const Property& wantedProperty) const {
  // Call this function several times to account for possible result updates throughout the autocompletion process
  PropertyList availableProperties = results.allContainedProperties();
  // Ensure that there is a mapping encoded for this property
  if (requiredPropertiesMap_.find(wantedProperty) != requiredPropertiesMap_.end()) {
    PropertyList requiredProperties = requiredPropertiesMap_.at(wantedProperty);
    return availableProperties.containsSubSet(requiredProperties);
  }
  return false;
}

void ResultsAutoCompleter::generateThermochemistry(Results& results, const AtomCollection& atomCollection) {
  // Get spin multiplicity
  LcaoUtils::ElectronicOccupation electronicOccupation = results.get<Property::ElectronicOccupation>();
  int spinMultiplicity = abs(electronicOccupation.numberAlphaElectrons() - electronicOccupation.numberBetaElectrons()) + 1;
  generateThermochemistry(results, atomCollection, spinMultiplicity);
}

void ResultsAutoCompleter::generateThermochemistry(Results& results, const AtomCollection& atomCollection, int spinMultiplicity) {
  double electronicEnergy = results.get<Property::Energy>();
  std::unique_ptr<ThermochemistryCalculator> thermochemistryCalculator =
      std::make_unique<ThermochemistryCalculator>(results.get<Property::Hessian>(), atomCollection.getElements(),
                                                  atomCollection.getPositions(), spinMultiplicity, electronicEnergy);
  thermochemistryCalculator->setTemperature(temperature_);
  // TODO: When available incorporate molassembler's symmetry detection here
  thermochemistryCalculator->setMolecularSymmetryNumber(sigma_);
  thermochemistryCalculator->setZPVEInclusion(zpveIncluded_);
  ThermochemicalComponentsContainer thermochemicalComponentsContainer = thermochemistryCalculator->calculate();
  results.set<Property::Thermochemistry>(std::move(thermochemicalComponentsContainer));
}

void ResultsAutoCompleter::generateDensityMatrix(Scine::Utils::Results& results) {
  DensityMatrix densityMatrix = LcaoUtils::DensityMatrixGenerator::generate(results.get<Property::ElectronicOccupation>(),
                                                                            results.get<Property::CoefficientMatrix>());
  results.set<Property::DensityMatrix>(std::move(densityMatrix));
}

void ResultsAutoCompleter::generateAtomicCharges(Results& results) {
  // Initialize atomic charges to vector of nAtoms zeros
  std::vector<double> atomicCharges(coreCharges_.size());
  LcaoUtils::calculateMullikenAtomicCharges(atomicCharges, coreCharges_, results.get<Property::DensityMatrix>(),
                                            results.get<Property::OverlapMatrix>(),
                                            results.get<Property::AOtoAtomMapping>());
  results.set<Property::AtomicCharges>(std::move(atomicCharges));
}

void ResultsAutoCompleter::generateBondOrderMatrix(Results& results) {
  BondOrderCollection bondOrderMatrix(coreCharges_.size());
  LcaoUtils::calculateBondOrderMatrix(bondOrderMatrix, results.get<Property::DensityMatrix>(),
                                      results.get<Property::OverlapMatrix>(), results.get<Property::AOtoAtomMapping>());
  results.set<Property::BondOrderMatrix>(std::move(bondOrderMatrix));
}

void ResultsAutoCompleter::generateProperties(Results& results, const AtomCollection& atomCollection) {
  // To ensure that also properties derivable from other derived properties are calculated (even if the calculation
  // routines are not called in the optimal order here) repeat until no new properties could be derived anymore
  bool calculatedNewProperty = true;
  while (calculatedNewProperty) {
    calculatedNewProperty = false;
    // Loop over properties and check whether they are wanted
    for (Property propertyToCompute : allProperties) {
      if (wantedProperties_.containsSubSet(propertyToCompute)) {
        // Ensure that property does not yet exist to avoid double effort
        PropertyList availableProperties = results.allContainedProperties();
        if (not availableProperties.containsSubSet(propertyToCompute)) {
          // Check whether property can be generated
          if (this->propertyGeneratable(results, propertyToCompute)) {
            calculatedNewProperty = true;
            // TODO Think about more elegant/easier to maintain set-up! Map to function pointers?
            // Compute property and store it in results
            if (propertyToCompute == Property::Thermochemistry) {
              this->generateThermochemistry(results, atomCollection);
            }
            else if (propertyToCompute == Property::DensityMatrix) {
              generateDensityMatrix(results);
            }
            else if (propertyToCompute == Property::AtomicCharges) {
              this->generateAtomicCharges(results);
            }
            else if (propertyToCompute == Property::BondOrderMatrix) {
              this->generateBondOrderMatrix(results);
            }
            else {
              throw std::logic_error("No calculation routine provided for Property " +
                                     std::to_string(static_cast<std::underlying_type<Property>::type>(propertyToCompute)));
            }
          } // If: property generatable
        }   // If not: property already in results
      }     // If: property wanted
    }       // For: over all properties
  }         // While: new property could be derived
}

} // namespace Utils
} // namespace Scine
