/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_GAUSSIANCALCULATOR_H
#define UTILS_GAUSSIANCALCULATOR_H

#include <Core/Interfaces/Calculator.h>
#include <Core/Interfaces/ObjectWithOrbitals.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/DataStructures/MolecularOrbitals.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Settings.h>
#include <Utils/Technical/CloneInterface.h>

namespace Scine {
namespace Utils {
namespace ExternalQC {

class GaussianCalculator final : public Utils::CloneInterface<GaussianCalculator, Core::Calculator>,
                                 public Core::ObjectWithOrbitals {
 public:
  static constexpr const char* model = "GAUSSIAN";
  /// @brief Default Constructor.
  GaussianCalculator();
  /// @brief Default Destructor.
  ~GaussianCalculator() final = default;
  /// @brief Copy Constructor.
  GaussianCalculator(const GaussianCalculator& rhs);
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
   * @brief Getter for the properties to calculate.
   */
  PropertyList getRequiredProperties() const override;
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
   * @brief Updates the molecular orbitals stored in the checkpoint file
   * @param mos The molecular orbitals
   * @note This requires a checkpoint file to be present in the calculation directory
   */
  void setOrbitals(const MolecularOrbitals& mos) override;

  /**
   * @brief Getter for the name of the Calculator.
   * @return Returns the name of the Calculator.
   */
  std::string name() const override;
  /**
   * @brief Accessor for the settings.
   * @return Utils::Settings& The settings.
   */
  Utils::Settings& settings() final;
  /**
   * @brief Constant accessor for the settings.
   * @return const Utils::Settings& The settings.
   */
  const Utils::Settings& settings() const final;
  /**
   * @brief Implements Core::StateHandableObject::getState().
   * Note: Not implemented, will throw an exception.
   * @return std::shared_ptr<Core::State> The current state
   */
  std::shared_ptr<Core::State> getState() const final;
  /**
   * @brief Implements Core::StateHandableObject::loadState().
   * Note: Not implemented, will throw an exception.
   * @param state The new state.
   */
  void loadState(std::shared_ptr<Core::State> state) final;
  /**
   * @brief Accessor for the saved instance of Utils::Results.
   * @return Utils::Results& The results of the previous calculation.
   */
  Utils::Results& results() final;
  /**
   * @brief Constant accessor for the Utils::Results.
   * @return const Utils::Results& The results of the previous calculation.
   */
  const Utils::Results& results() const final;
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
  /**
   * @brief Whether the calculator has no underlying Python code and can therefore
   * release the global interpreter lock in Python bindings
   */
  bool allowsPythonGILRelease() const override {
    return true;
  };

 private:
  /*
   * @brief Checks whether the given Gaussian binary is valid.
   */
  bool binaryIsValid();
  /*
   * @brief Implementation of a calculation.
   */
  const Utils::Results& calculateImpl(std::string description);
  /*
   * @brief Apply settings.
   */
  void applySettings();
  // The settings.
  std::unique_ptr<Utils::Settings> settings_;
  // The results.
  Utils::Results results_;
  // Base working directory
  std::string baseWorkingDirectory_;
  // Calculation directory
  std::string calculationDirectory_;
  // Settings as private members:
  std::string fileNameBase_;
  std::string gaussianExecutable_ = "";
  std::string gaussianDirectory_ = "";
  AtomCollection atoms_;
  Utils::PropertyList requiredProperties_;
  // Keeps track of whether the binary has been checked for validity yet
  bool binaryHasBeenChecked_ = false;
  const std::vector<std::string> availableSolvationModels_ =
      std::vector<std::string>{"cpcm", "pcm", "dipole", "ipcm", "scipcm", "smd"};
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_GAUSSIANCALCULATOR_H
