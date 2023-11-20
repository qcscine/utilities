/**
 * @file
 * @brief Default module that is always loaded by Core
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_MRCCDFTCALCULATOR_H
#define UTILS_MRCCDFTCALCULATOR_H

#include "Utils/ExternalQC/MRCC/MrccCalculator.h"

namespace Scine::Utils::ExternalQC {

/**
 * @class
 * @brief DFT calculator for MRCC.
 */
class MrccDFTCalculator final
  : public Scine::Utils::CloneInterface<MrccDFTCalculator, MrccCalculator, Scine::Core::Calculator> {
 public:
  MrccDFTCalculator() = default;
  ~MrccDFTCalculator() = default;
  static constexpr const char* model = "DFT";

 private:
  std::string getMethodFamily() const override final;
};

} // namespace Scine::Utils::ExternalQC

#endif // UTILS_MRCCDFTCALCULATOR_H
