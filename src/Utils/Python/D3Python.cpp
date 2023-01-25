/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Dftd3/Dftd3.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

using namespace Scine::Utils;

void init_and_calculate(Dftd3::Dftd3& self, const AtomCollection& atomCollection, double s6, double s8,
                        double dampingParam1, double dampingParam2, Dftd3::Damping damping = Dftd3::Damping::BJ) {
  self.initialize(atomCollection, s6, s8, dampingParam1, dampingParam2, damping);
  self.calculate(Derivative::First);
}

void init_d3(pybind11::module& m) {
  pybind11::enum_<Dftd3::Damping> damping(m, "Damping");
  damping.value("BJ", Dftd3::Damping::BJ);
  damping.value("Zero", Dftd3::Damping::Zero);

  pybind11::class_<Dftd3::Dftd3> disp(m, "D3Evaluator");
  disp.def(pybind11::init<>());
  disp.def("calculate", &init_and_calculate,
           "Arguments:"
           "atomCollection The atom collection (molecule) for which the D3 correction should be calculated."
           "s6 The s6 scaling parameter."
           "s8 The s8 scaling parameter."
           "dampingParam1 The first parameter of the damping function (a1 for BJ damping, sr for zero damping)."
           "dampingParam2 The second parameter of the damping function (a2 for BJ damping, a for zero damping)."
           "damping The damping function that should be used.",
           pybind11::arg("atomCollection"), pybind11::arg("s6"), pybind11::arg("s8"), pybind11::arg("dampingParam1"),
           pybind11::arg("dampingParam2"), pybind11::arg("damping") = Dftd3::Damping::BJ);
  disp.def("get_energy", &Dftd3::Dftd3::getEnergy, "Getter for the D3 energy correction.");
  disp.def("get_gradients", &Dftd3::Dftd3::getGradients, "Getter for the D3 nuclear gradient correction.");
}
