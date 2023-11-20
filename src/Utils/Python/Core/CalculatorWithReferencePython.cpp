/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
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

/**
 * Adapted from Calculator python bindings.
 */

namespace PythonHelperFunctions {
inline std::shared_ptr<Scine::Core::CalculatorWithReference> getCalculatorWithReference(std::string name, std::string program) {
  // Load module manager
  auto& manager = Scine::Core::ModuleManager::getInstance();

  // Transform to standard case: Program and NAME
  for (auto& x : program)
    x = std::tolower(x);
  program[0] = std::toupper(program[0]);
  for (auto& x : name)
    x = std::toupper(x);

  // Generate Calculator
  std::shared_ptr<Scine::Core::CalculatorWithReference> calc;
  try {
    calc = manager.get<Scine::Core::CalculatorWithReference>(name, program);
  }
  catch (...) {
    if (program.empty()) {
      std::cout << "No SCINE module providing '" << name << "' is currently loaded.\n";
      std::cout << "Please add the module to the SCINE_MODULE_PATH\n";
      std::cout << "or load the corresponding Python module in order for it to be accessible.\n";
      throw std::runtime_error("Failed to load method/program.");
    }
    else {
      std::cout << "No SCINE module named '" << program << "' providing '" << name << "' is currently loaded.\n";
      std::cout << "Please add the module to the SCINE_MODULE_PATH\n";
      std::cout << "or load the corresponding Python module in order for it to be accessible.\n";
      throw std::runtime_error("Failed to load method/program.");
    }
  }
  // Return CalculatorWithReference
  return calc;
}

inline bool hasResults(const Scine::Core::CalculatorWithReference& calculator) {
  try {
    calculator.results();
  }
  catch (...) {
    return false;
  }
  return true;
}

inline std::vector<std::string> listSettings(const Scine::Utils::Settings& settings) {
  std::vector<std::string> availableSettings;
  for (const auto& s : settings.getDescriptorCollection()) {
    availableSettings.push_back(s.first);
  }
  return availableSettings;
}

} // namespace PythonHelperFunctions

void init_calculator_with_reference(pybind11::module& m) {
  pybind11::class_<Scine::Core::CalculatorWithReference, std::shared_ptr<Scine::Core::CalculatorWithReference>> calculatorWithReference(
      m, "CalculatorWithReference");

  calculatorWithReference.doc() = "The CalculatorWithReference is the abstract base for classes running calculations "
                                  "based on a reference one, f.i. Excited-state calculations.";

  /* NOTE: This is absolutely necessary since Core doesn't STORE the constexpr
   * string anywhere, but this way there should be no missing symbols
   */
  constexpr const char* calcWithRefInterfaceStr = Scine::Core::CalculatorWithReference::interface;
  calculatorWithReference.attr("INTERFACE") = calcWithRefInterfaceStr;

  calculatorWithReference.def_property(
      "reference_calculator",
      pybind11::overload_cast<>(&Scine::Core::CalculatorWithReference::getReferenceCalculator, pybind11::const_),
      &Scine::Core::CalculatorWithReference::setReferenceCalculator, "The reference calculator.");

  calculatorWithReference.def_property(
      "settings", [](Scine::Core::CalculatorWithReference& calc) -> Scine::Utils::Settings& { return calc.settings(); },
      [](Scine::Core::CalculatorWithReference& calc, Scine::Utils::Settings settings) {
        calc.settings() = std::move(settings);
      },
      pybind11::return_value_policy::reference, "Settings of the calculator");

  calculatorWithReference.def_property(
      "log", [](Scine::Core::CalculatorWithReference& self) -> Scine::Core::Log& { return self.getLog(); },
      [](Scine::Core::CalculatorWithReference& self, Scine::Core::Log log) -> void { self.setLog(std::move(log)); },
      pybind11::return_value_policy::reference, "Logger of the calculator with reference.");

  calculatorWithReference.def("apply_settings", &Scine::Core::CalculatorWithReference::applySettings,
                              "Applies the settings given.");

  calculatorWithReference.def("reference_calculation", &Scine::Core::CalculatorWithReference::referenceCalculation,
                              "Execute the reference calculation");
  calculatorWithReference.def("calculate", &Scine::Core::CalculatorWithReference::calculate,
                              "Execute the calculation. Needs a reference_calculation call beforehand.");

  calculatorWithReference.def("get_results", pybind11::overload_cast<>(&Scine::Core::CalculatorWithReference::results),
                              "Get any strored results.");
  calculatorWithReference.def("has_results", &PythonHelperFunctions::hasResults, "Check if results are present.");

  calculatorWithReference.def("name", &Scine::Core::CalculatorWithReference::name, "Yields the name of the calculator");

  calculatorWithReference.def(
      "list_settings",
      [](const Scine::Core::CalculatorWithReference& self) -> std::vector<std::string> {
        return PythonHelperFunctions::listSettings(self.settings());
      },
      "Returns the list of settings available in this instance.");
}
