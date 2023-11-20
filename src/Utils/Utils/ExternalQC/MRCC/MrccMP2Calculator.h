/**
 * @file
 * @brief Default module that is always loaded by Core
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_MRCCMP2CALCULATOR_H
#define UTILS_MRCCMP2CALCULATOR_H

#include "Utils/ExternalQC/MRCC/MrccCalculator.h"

namespace Scine::Utils::ExternalQC {

/**
 * @class
 * @brief MP2 calculator for MRCC.
 */
class MrccMP2Calculator final
  : public Scine::Utils::CloneInterface<MrccMP2Calculator, MrccCalculator, Scine::Core::Calculator> {
 public:
  MrccMP2Calculator() = default;
  ~MrccMP2Calculator() = default;
  static constexpr const char* model = "MP2";

 private:
  std::string getMethodFamily() const override final;
};

} // namespace Scine::Utils::ExternalQC

#endif // UTILS_MRCCMP2CALCULATOR_H
