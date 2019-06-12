/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/Geometry/AtomCollection.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

using namespace Scine::Core;

void init_calculator(pybind11::module& m) {
  pybind11::class_<Calculator, std::shared_ptr<Calculator>> calculator(m, "Calculator");

  calculator.def("set_structure", &Calculator::setStructure, "Set the molecular structure to calculate");

  calculator.def_property("positions", &Calculator::getPositions, &Calculator::modifyPositions,
                          "Positions of the molecular structure");

  calculator.def("set_required_properties", &Calculator::setRequiredProperties, "Set the required properties");

  calculator.def("possible_properties", &Calculator::possibleProperties,
                 "Yields a list of properties that the calculator can calculate");

  calculator.def("calculate", &Calculator::calculate, pybind11::arg("dummy") = std::string{}, "Execute the calculation");

  calculator.def("name", &Calculator::name, "Yields the name of the calculator");

  // TODO missing: settings, statesHandler (these need to be bound for them to be useful)
}
