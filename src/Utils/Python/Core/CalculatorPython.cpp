/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Core/Interfaces/Calculator.h>
#include <Core/Interfaces/CalculatorWithReference.h>
#include <Core/ModuleManager.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/Geometry.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/IO/Yaml.h>
#include <Utils/Settings.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <yaml-cpp/yaml.h>
#include <iostream>

std::shared_ptr<Scine::Core::Calculator> getCalculator(std::string method_family, std::string program) {
  // Load module manager
  auto& manager = Scine::Core::ModuleManager::getInstance();

  for (auto& x : program)
    x = std::tolower(x);
  program[0] = std::toupper(program[0]);
  for (auto& x : method_family)
    x = std::toupper(x);

  // Generate Calculator
  std::shared_ptr<Scine::Core::Calculator> calc;
  try {
    calc = manager.get<Scine::Core::Calculator>(Scine::Core::Calculator::supports(method_family), program);
  }
  catch (...) {
    if (program.empty()) {
      std::cout << "No SCINE module providing '" << method_family << "' is currently loaded.\n";
      std::cout << "Please add the module to the SCINE_MODULE_PATH\n";
      std::cout << "or load the corresponding Python module in order for it to be accessible.\n";
      throw std::runtime_error("Failed to load method/program.");
    }
    else {
      std::cout << "No SCINE module named '" << program << "' providing '" << method_family << "' is currently loaded.\n";
      std::cout << "Please add the module to the SCINE_MODULE_PATH\n";
      std::cout << "or load the corresponding Python module in order for it to be accessible.\n";
      throw std::runtime_error("Failed to load method/program.");
    }
  }
  // Return Calculator
  return calc;
}

std::shared_ptr<Scine::Core::Calculator> loadSystem(std::string path, std::string method_family, pybind11::kwargs kwargs) {
  // Read all optional arguments and turn them into a node
  std::string program = "";
  YAML::Node settingsnode;
  for (auto item : kwargs) {
    std::string key = item.first.cast<std::string>();
    std::string value = pybind11::str(item.second).cast<std::string>();
    if (!key.compare("program")) {
      program = value;
    }
    else {
      settingsnode[key] = value;
    }
  }

  // Load molecule
  auto readResults = Scine::Utils::ChemicalFileHandler::read(path);
  // Generate Calculator
  auto calc = getCalculator(method_family, program);
  // Apply settings to Calculator
  nodeToSettings(calc->settings(), settingsnode);
  // Set initial structure
  calc->setStructure(readResults.first);
  return calc;
}

std::vector<std::string> getAvailableSettings(std::string method_family, std::string program) {
  // Get the correct Calculator
  auto calc = getCalculator(method_family, program);
  // Get the names of the available settings from the calculator's settings
  std::vector<std::string> availableSettings;
  Scine::Utils::UniversalSettings::DescriptorCollection settings = calc->settings().getDescriptorCollection();
  for (const auto& s : settings) {
    availableSettings.push_back(s.first);
  }
  return availableSettings;
}

Scine::Utils::PropertyList getPossiblePropertiesByStrings(std::string method_family, std::string program) {
  // Get the correct Calculator
  auto calc = getCalculator(method_family, program);
  return calc->possibleProperties();
}

void setRequiredProperties(Scine::Core::Calculator& calculator, std::vector<Scine::Utils::Property> list) {
  Scine::Utils::PropertyList properties;
  for (auto& p : list)
    properties.addProperty(p);
  calculator.setRequiredProperties(properties);
}

bool hasResults(const Scine::Core::Calculator& calculator) {
  try {
    calculator.results();
  }
  catch (...) {
    return false;
  }
  return true;
}

Scine::Utils::AtomCollection getStructure(const Scine::Core::Calculator& calc) {
  Scine::Utils::AtomCollection ret(*(calc.getStructure()));
  return ret;
}

void init_calculator(pybind11::module& m) {
  pybind11::class_<Scine::Core::Calculator, std::shared_ptr<Scine::Core::Calculator>> calculator(m, "Calculator");

  calculator.doc() = "The Calculator is the abstract base for classes running electronic structure calculations.";

  /* NOTE: This is absolutely necessary since Core doesn't STORE the constexpr
   * string anywhere, but this way there should be no missing symbols
   */
  constexpr const char* calcInterfaceStr = Scine::Core::Calculator::interface;
  calculator.attr("INTERFACE") = calcInterfaceStr;

  calculator.def_property("structure", &getStructure, &Scine::Core::Calculator::setStructure,
                          "The molecular structure to calculate");

  calculator.def_property("positions", &Scine::Core::Calculator::getPositions,
                          &Scine::Core::Calculator::modifyPositions, "Positions of the molecular structure");

  calculator.def_property(
      "settings", [](Scine::Core::Calculator& calc) -> Scine::Utils::Settings& { return calc.settings(); },
      [](Scine::Core::Calculator& calc, Scine::Utils::Settings settings) { calc.settings() = std::move(settings); },
      pybind11::return_value_policy::reference, "Settings of the calculator");

  calculator.def("set_required_properties", &Scine::Core::Calculator::setRequiredProperties, "Set the required properties");
  calculator.def("set_required_properties", &setRequiredProperties, "Set the required properties");

  calculator.def("get_possible_properties", &Scine::Core::Calculator::possibleProperties,
                 "Yields a list of properties that the calculator can calculate");

  calculator.def("calculate", &Scine::Core::Calculator::calculate, pybind11::arg("dummy") = std::string{},
                 "Execute the calculation");

  calculator.def("get_results", pybind11::overload_cast<>(&Scine::Core::Calculator::results), "Get any strored results.");
  calculator.def("has_results", &hasResults, "Check if results are present.");

  calculator.def("name", &Scine::Core::Calculator::name, "Yields the name of the calculator");

  // Some static helper functions
  m.def("load_system", &loadSystem, pybind11::arg("path"), pybind11::arg("method_family"),
        "Loads a single system (xyz-file) into a Calculator with the given method and optional settings. (Deprecated)");
  m.def("load_system_into_calculator", &loadSystem, pybind11::arg("path"), pybind11::arg("method_family"),
        "Loads a single system (xyz-file) into a Calculator with the given method and optional settings.");
  m.def("get_available_settings", &getAvailableSettings, pybind11::arg("method_family"), pybind11::arg("program"),
        "Lists the available settings of a Calculator with the given method and from the given program");
  m.def("get_possible_properties", &getPossiblePropertiesByStrings, pybind11::arg("method_family"), pybind11::arg("program"),
        "Lists the available properties of a Calculator with the given method and from the given program.");

  // TODO missing: statesHandler (these need to be bound for them to be useful)
}
