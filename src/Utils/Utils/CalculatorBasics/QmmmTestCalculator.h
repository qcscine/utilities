/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILSOS_QMMMTESTCALCULATOR_H
#define UTILSOS_QMMMTESTCALCULATOR_H

#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics/TestCalculator.h>

namespace Scine {
namespace Utils {

class QmmmTestCalculator : public Core::Calculator {
 public:
  QmmmTestCalculator();
  ~QmmmTestCalculator() override = default;
  void setStructure(const AtomCollection& structure) final;
  std::unique_ptr<Utils::AtomCollection> getStructure() const final;
  void modifyPositions(PositionCollection newPositions) final;
  const PositionCollection& getPositions() const final;
  void setRequiredProperties(const PropertyList& /* requiredProperties */) final;
  PropertyList getRequiredProperties() const final;
  PropertyList possibleProperties() const final;
  const Results& calculate(std::string dummy = "") final;
  std::string name() const final;
  std::shared_ptr<Core::State> getState() const final;
  void loadState(std::shared_ptr<Core::State> /* state */) final;
  const Settings& settings() const final;
  Settings& settings() final;
  Utils::Results& results() final;
  const Utils::Results& results() const final;
  bool supportsMethodFamily(const std::string& methodFamily) const final;
  bool allowsPythonGILRelease() const override;

 private:
  AtomCollection _structure;
  Results _results;
  std::unique_ptr<Settings> _settings;
  std::shared_ptr<TestCalculator> _underlyingCalculator;
  std::shared_ptr<Core::Calculator> cloneImpl() const override;
};

} // namespace Utils
} // namespace Scine

#endif // UTILSOS_QMMMTESTCALCULATOR_H
