/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_FOCK_SIMPLE_H
#define UTILS_FOCK_SIMPLE_H

#include <Utils/MethodEssentials/Methods/SCFModifier.h>
#include <Eigen/Core>
#include <vector>

namespace Scine {
namespace Utils {

class Fock_Simple : public SCFModifier {
 public:
  void onFockCalculated() override;

  void initialize() override;
  void setDamping(double damping) {
    damping_ = damping;
  }

  void addMatrices(const Eigen::MatrixXd& F);

 private:
  const Eigen::MatrixXd& extrapolate();
  bool initialized{false};
  double damping_{0.8};
  int nAOs_;
  int index_;
  std::vector<Eigen::MatrixXd> fockMatrices;
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_FOCK_SIMPLE_H
