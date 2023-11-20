/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_TESTCALCULATO_resultsH
#define UTILS_TESTCALCULATO_resultsH

#include <Core/Interfaces/Calculator.h>
#include <Core/Interfaces/WavefunctionOutputGenerator.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Technical/CloneInterface.h>

namespace Scine {
namespace Utils {

/**
 * @brief A test calculator using only an atom pairwise potential for energies
 *        and gradients.
 *
 * The exact energy function is
 *  E_tot = sum(i<j) E(r_ij),
 * with:
 *  E(r) = lj_depth(r_min/r)^12 - 2(r_min/r)^6 + g_height*exp(-((r-r_max)/g_width)^2)
 *
 * The distance dependent function is a Lennard-Jones potential with an
 * additional Gaussian added to generate repulsion of atoms that are not bonded.
 *
 * The pure Lennard-Jones potential defines a minimum at r_min while the
 * pure Gaussian generates a maximum at r_max interatomic distance.
 * However, due to the combination of the two, the actual values for the extrema
 * will be slightly different than the values entered into the function.
 *
 * All values are scaled and shifted depending on the two atoms in the given
 * pair. All shifts are based on the covalent radii of the two atoms. With
 * the reference being a C-H bond denoted with a scaling factor of 1.0.
 * For more details, see the implementation
 *
 * Analytical first derivatives and numerical second derivatives are available.
 *
 */
class TestCalculator : public CloneInterface<TestCalculator, Core::Calculator>, public Core::WavefunctionOutputGenerator {
 public:
  TestCalculator();
  ~TestCalculator() override = default;
  static constexpr const char* model = "TEST";
  bool supportsMethodFamily(const std::string& methodFamily) const override;
  void setStructure(const AtomCollection& structure) final;
  void modifyPositions(PositionCollection newPositions) final;
  const PositionCollection& getPositions() const final;
  void setRequiredProperties(const PropertyList& requiredProperties) final;
  PropertyList getRequiredProperties() const final;
  PropertyList possibleProperties() const final;
  const Results& calculate(std::string dummy = "") final;
  std::string name() const final;
  const Settings& settings() const final;
  Settings& settings() final;
  std::shared_ptr<Core::State> getState() const final;
  void loadState(std::shared_ptr<Core::State> state) final;
  Utils::Results& results() final;
  const Utils::Results& results() const final;
  std::unique_ptr<Utils::AtomCollection> getStructure() const final;
  void generateWavefunctionInformation(const std::string& out) final;
  void generateWavefunctionInformation(std::ostream& out) final;
  /**
   * @brief Set the precision of the calculator.
   * @param value The exponent to the expression 10^(-value).
   */
  void setPrecision(double value);
  /**
   * @brief Get the precision of the calculator.
   * @return The precision set for the calculator
   */
  double getPrecision();
  /**
   * @brief Whether the calculator has no underlying Python code and can therefore
   * release the global interpreter lock in Python bindings
   */
  bool allowsPythonGILRelease() const override {
    return true;
  };

 private:
  PropertyList _requiredProperties{};
  AtomCollection _structure;
  Results _results;
  std::shared_ptr<Settings> _settings;
  /**
   * @brief The default precision of the calculator.
   */
  double _precision = 13.0;
  /**
   * @brief Truncate a given value to the set precision.
   * @param value The value to be truncated
   * @return The truncated value.
   */
  double truncateOff(double value);
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_TESTCALCULATO_resultsH
