/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Math/IterativeDiagonalizer/SpinAdaptedEigenContainer.h>
#include <Utils/TimeDependent/TransitionDipoleCalculator.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;

void init_spin_adapted_electronic_transition_result(pybind11::module& m) {
  using RowMatrixXd = Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::RowMajor>;
  m.def("transition_dipole_to_oscillator_strength", [](const RowMatrixXd& m, const Eigen::VectorXd& eigenvalues) {
    return TransitionDipoleCalculator::transitionDipoleMomentToOscillatorStrength(Eigen::Matrix3Xd(m), eigenvalues);
  });

  pybind11::class_<EigenContainer, std::shared_ptr<EigenContainer>> eigenContainer(m, "EigenContainer");
  eigenContainer.def_readonly("eigenvalues", &EigenContainer::eigenValues);
  eigenContainer.def_readonly("eigenvectors", &EigenContainer::eigenVectors);

  pybind11::class_<ElectronicTransitionResult, std::shared_ptr<ElectronicTransitionResult>> electronicTransitionResult(
      m, "ElectronicTransitionResult");
  electronicTransitionResult.def_readonly("eigenstates", &ElectronicTransitionResult::eigenStates);
  electronicTransitionResult.def_readonly("transition_dipoles", &ElectronicTransitionResult::transitionDipoles);

  pybind11::class_<SpinAdaptedElectronicTransitionResult> spinAdaptedElectronicTransitionResult(
      m, "SpinAdaptedElectronicTransitionResult");
  spinAdaptedElectronicTransitionResult.def_readonly("singlet", &SpinAdaptedElectronicTransitionResult::singlet);
  spinAdaptedElectronicTransitionResult.def_readonly("triplet", &SpinAdaptedElectronicTransitionResult::triplet);
  spinAdaptedElectronicTransitionResult.def_readonly("unrestricted", &SpinAdaptedElectronicTransitionResult::unrestricted);
  spinAdaptedElectronicTransitionResult.def_readonly("mo_labels", &SpinAdaptedElectronicTransitionResult::transitionLabels);
}
