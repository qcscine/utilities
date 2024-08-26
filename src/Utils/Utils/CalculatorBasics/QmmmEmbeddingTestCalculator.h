/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILSOS_QMMMEMBEDDINGTESTCALCULATOR_H
#define UTILSOS_QMMMEMBEDDINGTESTCALCULATOR_H

#include "Core/Interfaces/Calculator.h"
#include "Core/Interfaces/EmbeddingCalculator.h"
#include "Utils/CalculatorBasics/TestCalculator.h"

namespace Scine {
namespace Utils {

class QmmmEmbeddingTestCalculator
  : public CloneInterface<QmmmEmbeddingTestCalculator, Core::EmbeddingCalculator, Core::Calculator> {
 public:
  QmmmEmbeddingTestCalculator();
  QmmmEmbeddingTestCalculator(const QmmmEmbeddingTestCalculator& rhs);
  ~QmmmEmbeddingTestCalculator() override = default;
  void setUnderlyingCalculators(std::vector<std::shared_ptr<Calculator>> underlyingCalculators) override;
  std::vector<std::shared_ptr<Core::Calculator>> getUnderlyingCalculators() const override;
  void addUnderlyingSettings() override;
  void setStructure(const AtomCollection& structure) override;
  std::unique_ptr<Utils::AtomCollection> getStructure() const override;
  void modifyPositions(PositionCollection newPositions) override;
  const PositionCollection& getPositions() const override;
  void setRequiredProperties(const PropertyList& requiredProperties) override;
  PropertyList getRequiredProperties() const override;
  PropertyList possibleProperties() const override;
  const Results& calculate(std::string dummy = "") override;
  std::string name() const override;
  std::shared_ptr<Core::State> getState() const final;
  void loadState(std::shared_ptr<Core::State> /* state */) final;
  const Settings& settings() const override;
  Settings& settings() override;
  Utils::Results& results() override;
  const Utils::Results& results() const override;
  bool supportsMethodFamily(const std::string& methodFamily) const override;
  bool allowsPythonGILRelease() const override;

 private:
  void setStructureImpl(const Utils::AtomCollection& structure);
  void setUnderlyingCalculatorsImpl(std::vector<std::shared_ptr<Core::Calculator>> underlyingCalculators);
  const Utils::Results& calculateImpl(std::string description);

  AtomCollection structure_;
  AtomCollection qmRegion_;
  Results results_;
  std::unique_ptr<Settings> settings_;
  Utils::PropertyList requiredProperties_;
  std::shared_ptr<Core::Calculator> qmCalculator_;
  std::shared_ptr<Core::Calculator> mmCalculator_;
};

} // namespace Utils
} // namespace Scine

#endif // UTILSOS_QMMMEMBEDDINGTESTCALCULATOR_H
