/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_OPEN_SOURCE_AUTOCOMPLETION_H
#define UTILS_OPEN_SOURCE_AUTOCOMPLETION_H

#include <Utils/CalculatorBasics/Results.h>
#include <Utils/Geometry/AtomCollection.h>
#include <map>

namespace Scine {
namespace Utils {

/**
 * @class ResultsAutoCompleter ResultsAutoCompleter.h
 * @brief Provides functionality to autocomplete the results of a calculation.
 */
class ResultsAutoCompleter {
 public:
  /**
   * @brief Generates an instance of ResultsAutoCompleter.
   *
   * The Density Matrix, Thermochemistry, Atomic Charges and Bond Orders will be set as target properties. This
   * selection can be changed e.g. via setWantedProperties(). Use setTemperature(), setMolecularSymmetryNumber() and
   * setZPVEInclusion() to adapt the settings used for calculating thermochemical properties. The core charges used
   * during the computation are set equal to the nuclear charges as a default and can be adapted via setCoreCharges().
   */
  ResultsAutoCompleter(AtomCollection& atomCollection);

  /**
   * @brief Adds wantedProperty to the list of properties to be autocompleted.
   *
   * Note that this does not guarantee that the
   * property can be generated from the available results (see propertyGeneratable()).
   */
  void addOneWantedProperty(Property wantedProperty);

  /**
   * @brief Sets the list of properties wanted to be autocompleted to PropertyList wantedProperties.
   *
   * Note that this does not guarantee that all of these properties can be generated from the results (see
   * propertyGeneratable()).
   */
  void setWantedProperties(PropertyList wantedProperties);

  /**
   * @brief Sets all properties encoded in PropertyList as targets for autocompletion.
   *
   * Note that this does not guarantee that all of them can be generated (see propertyGeneratable()).
   */
  void setAllPropertiesAsWanted();

  /**
   * @brief Checks whether property wantedProperty can be autocompleted from the available results
   *
   * The following properties can be derived if the listed results are present:
   *
   * Derived Property | Necessary Result Properties
   * -----------------|---------------------------------------------------------------------
   * %DensityMatrix   | SuccessfulCalculation, ElectronicOccupation, CoefficientMatrix
   * Thermochemistry  | SuccessfulCalculation, Hessian, ElectronicOccupation
   * AtomicCharges    | SuccessfulCalculation, %DensityMatrix, OverlapMatrix, AOtoAtomMapping
   * BondOrderMatrix  | SuccessfulCalculation, %DensityMatrix, OverlapMatrix, AOtoAtomMapping
   *
   * Properties for which no dependencies are given are always considered to be not generatable.
   */

  bool propertyGeneratable(const Results& results, const Property& wantedProperty) const;

  /**
   * @brief Computes properties from wantedPropertyList from results if possible according to propertyGeneratable()
   *
   * For all properties in the list of wanted properties it is checked, whether they can be computed from the properties
   * already present in the results. If yes, the property is calculated and added to the results.
   * Since the results are updated during the run, initially ungeneratable properties may still be computed during the
   * run if that requires only the prior computation of other derived and generatable properties.
   */
  void generateProperties(Results& results, const AtomCollection& atomCollection);

  /**
   * @brief Sets the temperature employed in the ThermochemistryCalculator class.
   *
   * The default is 298.15 K.
   */
  void setTemperature(double temperature);

  /**
   * @brief Sets the molecular symmetry number employed in the ThermochemistryCalculator.
   *
   * The default is 1. Employing a
   * wrong symmetry number results into an incorrect rotational entropy contribution.
   */
  void setMolecularSymmetryNumber(int sigma);

  /**
   * @brief Sets whether the energy already included the zero point vibrational energy.
   *
   * The default is that it is not included.
   */
  void setZPVEInclusion(ZPVEInclusion zpveInclusion);

  /**
   * @brief Sets the core charges to the nuclear charges.
   */
  void setCoreChargesToNuclearCharges(AtomCollection& atomCollection);

  /**
   * @brief Sets the core Charges employed in the calculation of AtomicCharges.
   *
   * They default to be the nuclear charges,
   * but have to be adapted dependent on the employed method, e.g., when effective core potentials are used.
   */
  void setCoreCharges(std::vector<double> coreCharges);

 private:
  std::map<Property, PropertyList> requiredPropertiesMap_;
  /** \brief  List of Properties to be autocompleted if possible */
  PropertyList wantedProperties_;
  /** \brief Temperature employed in ThermochemistryCalculator */

  double temperature_{298.15};
  int sigma_{1};
  ZPVEInclusion zpveIncluded_{ZPVEInclusion::notIncluded};
  std::vector<double> coreCharges_;

  /**
   * @brief Sets requiredPropertiesMap_ to hardcoded values.
   */
  void fillRequiredPropertiesMap();

  /**
   * @brief Sets the default choice for properties to be generated from the results if possible.
   */
  void setDefaultWantedProperties();

  /**
   * @brief Calculates Thermochemistry from results and stores it there.
   *
   * The default symmetry factor (1) is assumed, potentially leading to wrong rotational contributions for symmetric
   * molecules.
   * The default temperature of 298.5 K is assumed.
   *
   */
  void generateThermochemistry(Results& results, const AtomCollection& atomCollection);

  /**
   * @brief Calculates DensityMatrix from results and stores it there.
   */
  static void generateDensityMatrix(Results& results);

  /**
   * @brief Calculates AtomicCharges (Mulliken charges) from results and stores them there.
   */
  void generateAtomicCharges(Results& results);

  /**
   * @brief Calculates BondOrderMatrix (Mayer bond orders) from results and stores them there.
   */
  void generateBondOrderMatrix(Results& results);
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_OPEN_SOURCE_AUTOCOMPLETION_H
