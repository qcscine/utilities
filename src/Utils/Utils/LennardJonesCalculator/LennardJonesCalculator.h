/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_LENNARDJONESCALCULATOR_H
#define UTILS_LENNARDJONESCALCULATOR_H

#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/DataStructures/PeriodicBoundaries.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Settings.h>
#include <Utils/Technical/CloneInterface.h>

namespace Scine {
namespace Utils {

/**
 * @brief A Lennard-Jones calculator with periodic boundary conditions in a cubic box.
 *
 * The energy function is <BR>
 * \f$ E_{tot} = \sum_{\text{i}<\text{j}} E(r_\text{ij})\f$ <BR>
 * \f$ E(r) = \begin{cases} 4 \epsilon [\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6] & r <
 * r_\text{cut-off}\\ 0 & r \geq r_\text{cut-off} \end{cases} \f$
 *
 * @note This Lennard-Jones calculator has a hard cut-off without a shift or a tail correction being applied.
 */
class LennardJonesCalculator : public CloneInterface<LennardJonesCalculator, Core::Calculator> {
 public:
  LennardJonesCalculator();
  ~LennardJonesCalculator() override = default;
  static constexpr const char* model = "LENNARDJONES";
  bool supportsMethodFamily(const std::string& methodFamily) const override;
  void setStructure(const AtomCollection& structure) final;
  void modifyPositions(PositionCollection newPositions) final;
  const PositionCollection& getPositions() const final;
  void setRequiredProperties(const PropertyList& requiredProperties) final;
  PropertyList getRequiredProperties() const final;
  PropertyList possibleProperties() const final;
  const Results& calculate(std::string description) final;
  std::string name() const final;
  const Settings& settings() const final;
  Settings& settings() final;
  std::shared_ptr<Core::State> getState() const final;
  void loadState(std::shared_ptr<Core::State> state) final;
  Utils::Results& results() final;
  const Utils::Results& results() const final;
  std::unique_ptr<Utils::AtomCollection> getStructure() const final;

 private:
  /*
   * @brief Apply settings.
   */
  void applySettings();
  /*
   * @brief Calculate vector between p1 and p2 taking into account optional PBCs
   *
   * The positions have to be moved into the cell before.
   */
  Displacement calculateDistanceVector(const Position& p1, Position& p2) const;
  PropertyList requiredProperties_{};
  AtomCollection structure_;
  Results results_;
  std::shared_ptr<Settings> settings_;
  double sigma_;   // sigma in bohr
  double epsilon_; // epsilon in Hartree
  double cutoff_;  // cut-off radius in bohr; Has to be smaller than half the box size
  std::shared_ptr<PeriodicBoundaries> pbc_ = nullptr; // The periodic boundary conditions
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_LENNARDJONESCALCULATOR_H
