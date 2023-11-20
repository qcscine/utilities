/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Core/Interfaces/Calculator.h>
#include <Core/Interfaces/CalculatorWithReference.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/MolecularDynamics/MolecularDynamics.h>
#include <Utils/Settings.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine;
using namespace Utils;

void init_molecular_dynamics(pybind11::module& m) {
  pybind11::class_<MolecularDynamics> molecular_dynamics(m, "MolecularDynamics");

  molecular_dynamics.def(pybind11::init<Core::Calculator&>(), pybind11::arg("calculator"),
                         R"delim(
                           Initialize the MolecularDynamics object.

                           :param calculator: The calculator with which the MD simulation is performed.
                         )delim");

  molecular_dynamics.def(pybind11::init<Core::CalculatorWithReference&>(), pybind11::arg("calculatorWithReference"),
                         R"delim(
                          Initialize the MolecularDynamics object.

                          :param calculatorWithReference: The calculator with reference calculator with which the MD simulation is performed.
                        )delim");

  molecular_dynamics.def("perform_md_simulation", &MolecularDynamics::performMDSimulation, pybind11::arg("structure"),
                         pybind11::arg("logger"), R"delim(
                           Perform a molecular dynamics simulation.

                           :param structure: The initial molecular structure for the simulation.
                           :param logger: The logger.
                         )delim");

  molecular_dynamics.def("set_external_stop", &MolecularDynamics::setExternalStop, pybind11::arg("stop_function"),
                         R"delim(
                         Sets an external stop function for the MD. This function should return a boolean (whether to
                         terminate the MD simulation), and the arguments of the function are:
                         PositionCollection, Results, int (step number).

                         :param stop function: The stop function.
                         )delim");

  molecular_dynamics.def("set_bias_potential", &MolecularDynamics::setBiasPotential, pybind11::arg("bias_potential"),
                         R"delim(
                         Sets an external bias potential for the MD. This function should return
                         two objects (as a tuple), an integer (bias energy) and a numpy array (bias gradients).
                         The arguments of the function are: PositionCollection, Results, int (step number).

                         :param bias potential: The bias potential function.
                         )delim");

  molecular_dynamics.def_property(
      "settings", [](MolecularDynamics& md) -> Scine::Utils::Settings& { return md.settings(); },
      [](MolecularDynamics& md, Scine::Utils::Settings settings) { md.settings() = std::move(settings); },
      pybind11::return_value_policy::reference, "Settings of the molecular dynamics simulation.");

  molecular_dynamics.def("get_molecular_trajectory", &MolecularDynamics::getMolecularTrajectory,
                         "Returns the molecular trajectory of the molecular dynamics simulation.");

  molecular_dynamics.def("get_velocities", &MolecularDynamics::getVelocities, R"delim(
    Gets velocities corresponding to the structures in the molecular trajectory of the MD simulation.

    :return: The velocities.
    )delim");

  molecular_dynamics.def("get_final_velocities", &MolecularDynamics::getFinalVelocities, R"delim(
    Gets the last velocities encountered during the MD simulation.

    These velocities can be different from the last element of the velocities obtainable via getVelocities() if the
    record frequency is not one.

    :return: The final velocities.
    )delim");

  molecular_dynamics.def("get_temperatures", &MolecularDynamics::getTemperatures, R"delim(
    Gets temperatures corresponding to the structures in the molecular trajectory of the MD simulation.

    :return: The temperatures in Kelvin.
    )delim");

  molecular_dynamics.def("set_initial_velocities", &MolecularDynamics::setInitialVelocities,
                         pybind11::arg("initial_velocities"),
                         R"delim(
                         Explicitly sets the initial velocities to be used for the MD simulation

                         :param initial_velocities: The initial velocities.
                         )delim");
}
