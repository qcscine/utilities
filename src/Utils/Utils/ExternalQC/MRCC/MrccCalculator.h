/**
 * @file
 * @brief Default module that is always loaded by Core
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILSOS_MRCCCALCULATOR_H
#define UTILSOS_MRCCCALCULATOR_H

#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Settings.h>
#include <Utils/Technical/CloneInterface.h>
#include <memory>
#include <string>

namespace Scine {
namespace Utils {
namespace ExternalQC {

/**
 * @class
 * @brief Base class for all MRCC calculators.
 */
class MrccCalculator : public CloneInterface<Scine::Utils::Abstract<MrccCalculator>, Core::Calculator> {
 public:
  static constexpr const char* program = "MRCC";
  static constexpr const char* binaryEnvVariable = "MRCC_BINARY_PATH";
  /// @brief Default Constructor.
  MrccCalculator();
  /// @brief Default Destructor.
  ~MrccCalculator() = default;
  /// @brief Copy Constructor.
  MrccCalculator(const MrccCalculator& rhs);
  /**
   * @brief Changes the molecular structure to calculate.
   * @param structure A new AtomCollection to save.
   */
  void setStructure(const AtomCollection& structure) override;
  /**
   * @brief Gets the molecular structure as a std::unique_ptr<AtomCollection>.
   * @return std::unique_ptr<AtomCollection>
   */
  std::unique_ptr<AtomCollection> getStructure() const override;
  /**
   * @brief Sets the properties to calculate.
   * @param requiredProperties A PropertyList, a sequence of bits that represent the
   *        properties that must be calculated.
   */
  void setRequiredProperties(const PropertyList& requiredProperties) override;
  /**
   * @brief Getter for the properties to calculate.
   */
  PropertyList getRequiredProperties() const override;
  /**
   * @brief Returns the list of the possible properties to calculate analytically.
   * By some method analytical hessian calculation is not possible. In this case the
   * hessian calculation is done seminumerically.
   */
  PropertyList possibleProperties() const override;
  /**
   * @brief The main function running calculations (dummy).
   *
   * @param description   The calculation description.
   * @return Result Return the result of the calculation. The object contains the
   *                       properties that were given as requirement by the
   *                       Calculator::setRequiredProperties function.
   */
  const Results& calculate(std::string description) override;
  /**
   * @brief Getter for the name of the Calculator.
   * @return Returns the name of the Calculator.
   */
  std::string name() const override;
  /**
   * @brief Accessor for the settings.
   * @return Settings& The settings.
   */
  Settings& settings() override;
  /**
   * @brief Accessor for the settings.
   * @return Settings& The settings.
   */
  const Settings& settings() const override;
  /**
   * @brief Accessor for the saved instance of Results.
   * @return Results& The results of the previous calculation.
   */
  Results& results() override;
  /**
   * @brief Constant accessor for the Results.
   * @return const Results& The results of the previous calculation.
   */
  const Results& results() const override;
  /**
   * @brief Whether the calculator supports a method family
   * @param methodFamily identifier for the method family
   * @return whether the calculator supports a method family
   */
  bool supportsMethodFamily(const std::string& methodFamily) const override;
  /**
   * @brief Allows to modify the positions of the underlying AtomCollection
   * @param newPositions the new positions to be assigned to the underlying AtomCollection
   */
  void modifyPositions(PositionCollection newPositions) override;
  /**
   * @brief Getter for the coordinates of the underlying AtomCollection
   */
  const PositionCollection& getPositions() const override;
  /**
   * @brief Implements Core::StateHandableObject::getState().
   * @return std::shared_ptr<Core::State> The current state
   */
  std::shared_ptr<Core::State> getState() const final;
  /**
   * @brief Implements Core::StateHandableObject::loadState().
   * @param state The new state.
   */
  void loadState(std::shared_ptr<Core::State> state) final;
  /**
   * @brief Getter for the calculation directory.
   */
  std::string getCalculationDirectory() const;
  /**
   * @breif Getter for the binary directory.
   */
  std::string getBinaryDirectory() const;
  /**
   * @brief Whether the calculator has no underlying Python code and can therefore
   * release the global interpreter lock in Python bindings
   */
  bool allowsPythonGILRelease() const override {
    return true;
  };

 private:
  void applySettings();
  const Results& calculateImpl(std::string description);
  /**
   * @brief Must be implemented in all realizations of this calculator and return the method family, e.g., DFT or HF.
   * @return The method family.
   */
  virtual std::string getMethodFamily() const = 0;

  const std::string name_ = "MRCC";
  std::string calculationDirectory_;
  std::string baseWorkingDirectory_;
  std::string binaryDirectory_;
  PropertyList requiredProperties_ = Utils::Property::Energy;
  std::unique_ptr<Settings> settings_;
  Results results_;
  const std::vector<std::string> availableSolvationModels_ = std::vector<std::string>{"iefpcm"};
  AtomCollection atoms_;
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILSOS_MRCCCALCULATOR_H
