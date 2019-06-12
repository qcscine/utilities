/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_ORCACALCULATOR_H
#define UTILS_ORCACALCULATOR_H

#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/CalculatorBasics/StatesHandler.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Settings.h>
#include <Utils/Technical/CloneInterface.h>

namespace Scine {
namespace Utils {
namespace ExternalQC {

class OrcaCalculator : public Utils::CloneInterface<OrcaCalculator, Core::Calculator> {
 public:
  static constexpr const char* model = "ORCA";
  /// @brief Default Constructor.
  OrcaCalculator();
  /// @brief Default Destructor.
  ~OrcaCalculator() final = default;
  /// @brief Copy Constructor.
  OrcaCalculator(const OrcaCalculator& rhs);
  /**
   * @brief Changes the molecular structure to calculate.
   * @param structure A new Utils::AtomCollection to save.
   */
  void setStructure(const Utils::AtomCollection& structure) override;
  /**
   * @brief Gets the molecular structure as a std::unique_ptr<Utils::AtomCollection>.
   * @return std::unique_ptr<Utils::AtomCollection>
   */
  std::unique_ptr<Utils::AtomCollection> getStructure() const override;
  /**
   * @brief Allows to modify the positions of the underlying Utils::AtomCollection
   * @param newPositions the new positions to be assigned to the underlying Utils::AtomCollection
   */
  void modifyPositions(Utils::PositionCollection newPositions) override;
  /**
   * @brief Getter for the coordinates of the underlying Utils::AtomCollection
   */
  const Utils::PositionCollection& getPositions() const override;
  /**
   * @brief Sets the properties to calculate.
   * @param requiredProperties A Utils::PropertyList, a sequence of bits that represent the
   *        properties that must be calculated.
   */
  void setRequiredProperties(const Utils::PropertyList& requiredProperties) override;
  /**
   * @brief Returns the list of the possible properties to calculate analytically.
   * By some method analytical hessian calculation is not possible. In this case the
   * hessian calculation is done seminumerically.
   */
  Utils::PropertyList possibleProperties() const override;
  /**
   * @brief The main function running calculations (dummy).
   *
   * @param description   The calculation description.
   * @return Utils::Result Return the result of the calculation. The object contains the
   *                       properties that were given as requirement by the
   *                       Calculator::setRequiredProperties function.
   */
  const Utils::Results& calculate(std::string description) override;
  /**
   * @brief Getter for the name of the Calculator.
   * @return Returns the name of the Calculator.
   */
  std::string name() const override;
  /**
   * @brief Accessor for the settings.
   * @return Utils::Settings& The settings.
   */
  Utils::Settings& settings() override;
  /**
   * @brief Constant accessor for the settings.
   * @return const Utils::Settings& The settings.
   */
  const Utils::Settings& settings() const override;
  /**
   * @brief Accessor for the StatesHandler.
   * @return Utils::StatesHandler& The StatesHandler.
   */
  Utils::StatesHandler& statesHandler() override;
  /**
   * @brief Constant accessor for the StatesHandler.
   * @return const Utils::StatesHandler& The StatesHandler.
   */
  const Utils::StatesHandler& statesHandler() const override;
  /**
   * @brief Accessor for the saved instance of Utils::Results.
   * @return Utils::Results& The results of the previous calculation.
   */
  Utils::Results& results() override;
  /**
   * @brief Constant accessor for the Utils::Results.
   * @return const Utils::Results& The results of the previous calculation.
   */
  const Utils::Results& results() const override;
  /**
   * @brief Getter for the file name base string.
   */
  std::string getFileNameBase();
  /**
   * @brief Getter for the calculation directory.
   */
  std::string getCalculationDirectory();

 private:
  /*
   * @brief Implementation of a calculation.
   */
  const Utils::Results& calculateImpl(std::string description);
  /*
   * @brief Apply settings.
   */
  void applySettings();
  /*
   * @brief Creates a random name for a calculation directory and creates this directory.
   */
  std::string createNameForCalculationDirectory();
  // The settings.
  std::unique_ptr<Utils::Settings> settings_;
  // The states handler.
  std::unique_ptr<Utils::StatesHandler> statesHandler_;
  // The results.
  Utils::Results results_;
  // Base working directory
  std::string baseWorkingDirectory_;
  // Calculation directory
  std::string calculationDirectory_;
  // Settings as private members:
  std::string fileNameBase_;
  std::string orcaExecutable_;
  AtomCollection atoms_;
  Utils::PropertyList requiredProperties_;
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_ORCACALCULATOR_H
