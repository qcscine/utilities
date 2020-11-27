/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_ORCACALCULATOR_H
#define UTILS_ORCACALCULATOR_H

#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Settings.h>
#include <Utils/Technical/CloneInterface.h>

namespace Scine {
namespace Utils {
// class ResultsAutoCompleter;

namespace ExternalQC {

class OrcaCalculator final : public CloneInterface<OrcaCalculator, Core::Calculator> {
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
   * @param structure A new AtomCollection to save.
   */
  void setStructure(const AtomCollection& structure) override;
  /**
   * @brief Gets the molecular structure as a std::unique_ptr<AtomCollection>.
   * @return std::unique_ptr<AtomCollection>
   */
  std::unique_ptr<AtomCollection> getStructure() const override;
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
   * @brief Constant accessor for the settings.
   * @return const Settings& The settings.
   */
  const Settings& settings() const override;
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
   * @brief Getter for the file name base string.
   */
  std::string getFileNameBase() const;
  /**
   * @brief Getter for the calculation directory.
   */
  std::string getCalculationDirectory() const;
  /**
   * @brief Whether the calculator supports a method family
   * @param methodFamily identifier for the method family
   * @return whether the calculator supports a method family
   */
  bool supportsMethodFamily(const std::string& methodFamily) const override;

 private:
  /*
   * @brief Checks whether the given Orca binary is valid.
   */
  bool binaryIsValid();
  /*
   * @brief Implementation of a calculation.
   */
  const Results& calculateImpl(std::string description);
  /*
   * @brief Apply settings.
   */
  void applySettings();
  /*
   * @brief Copies .gbw files.
   * @param from Base name of the .gbw file (only 'base' of path/to/base.gbw)
   * @param to New base name of the .gbw file (only 'base' of path/to/base.gbw)
   */
  void copyBackupFile(const std::string& from, const std::string& to) const;
  /*
   * @brief Creates a random name for a calculation directory and creates this directory.
   */
  std::string createNameForCalculationDirectory();
  /*
   * @brief Deletes all .tmp files in the calculation directory.
   */
  void deleteTemporaryFiles();
  // The settings.
  std::unique_ptr<Settings> settings_;
  // The results.
  Results results_;
  //  std::unique_ptr<ResultsAutoCompleter> autoCompleter_;
  // Base working directory
  std::string baseWorkingDirectory_;
  // Calculation directory
  std::string calculationDirectory_;
  // Settings as private members:
  std::string fileNameBase_;
  std::string orcaExecutable_ = "";
  AtomCollection atoms_;
  PropertyList requiredProperties_;
  // Keeps track of whether the binary has been checked for validity yet
  bool binaryHasBeenChecked_ = false;
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_ORCACALCULATOR_H
