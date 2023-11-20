/**
 * @file
 * @brief Default module that is always loaded by Core
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_MRCCHFCALCULATOR_H
#define UTILS_MRCCHFCALCULATOR_H

#include "Utils/ExternalQC/MRCC/MrccCalculator.h"

namespace Scine::Utils::ExternalQC {

/**
 * @class
 * @brief Hartree--Fock calculator for MRCC.
 */
class MrccHFCalculator final : public Scine::Utils::CloneInterface<MrccHFCalculator, MrccCalculator, Scine::Core::Calculator> {
 public:
  MrccHFCalculator() = default;
  ~MrccHFCalculator() = default;
  static constexpr const char* model = "HF";

 private:
  std::string getMethodFamily() const override final;
};

} // namespace Scine::Utils::ExternalQC

#endif // UTILS_MRCCHFCALCULATOR_H
