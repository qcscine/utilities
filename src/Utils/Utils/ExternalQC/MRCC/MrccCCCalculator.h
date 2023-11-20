/**
 * @file
 * @brief Default module that is always loaded by Core
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_MRCCCCCALCULATOR_H
#define UTILS_MRCCCCCALCULATOR_H

#include "Utils/ExternalQC/MRCC/MrccCalculator.h"

namespace Scine::Utils::ExternalQC {

/**
 * @class
 * @brief Coupled cluster calculator for MRCC.
 */
class MrccCCCalculator final : public Scine::Utils::CloneInterface<MrccCCCalculator, MrccCalculator, Scine::Core::Calculator> {
 public:
  MrccCCCalculator() = default;
  ~MrccCCCalculator() = default;
  static constexpr const char* model = "CC";

 private:
  std::string getMethodFamily() const override final;
};

} // namespace Scine::Utils::ExternalQC

#endif // UTILS_MRCCCCCALCULATOR_H
