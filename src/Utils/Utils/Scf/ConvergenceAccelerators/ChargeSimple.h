/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_CHARGE_SIMPLE_H
#define UTILS_CHARGE_SIMPLE_H

#include <Utils/Scf/MethodInterfaces/ScfModifier.h>
#include <vector>

namespace Scine {
namespace Utils {

class ChargeSimple : public ScfModifier {
 public:
  void initialize() override;

  void setDamping(double damping) {
    damping_ = damping;
  }

  void addVector(const std::vector<double>& q);
  const std::vector<double>& extrapolate();

  void onIterationStart() override;

 private:
  double damping_{0.8};
  bool initialized{false};
  int nAtoms_;
  int index_;
  std::vector<std::vector<double>> qVectors;
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_CHARGE_SIMPLE_H
