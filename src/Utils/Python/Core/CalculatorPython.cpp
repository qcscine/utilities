/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Core/Interfaces/Calculator.h>
#include <Core/Interfaces/EmbeddingCalculator.h>
#include <Core/ModuleManager.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/Geometry.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/ValueCollection.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <algorithm>
#include <iostream>

// Defined in ValueCollectionPython.cpp
extern void update(Scine::Utils::UniversalSettings::ValueCollection& coll, const pybind11::dict& dictionary, bool preserveTypes);

void inputPreparation(std::string& method_family, std::string& program) {
  if (!program.empty()) {
    std::transform(std::begin(program) + 1, std::end(program), std::begin(program) + 1,
                   [](unsigned char x) { return std::tolower(x); });

    program[0] = std::toupper(program[0]);
  }

  std::transform(std::begin(method_family), std::end(method_family), std::begin(method_family),
                 [](unsigned char x) { return std::toupper(x); });
}

bool hasCalculator(std::string method_family, std::string program) {
  // Load module manager
  auto& manager = Scine::Core::ModuleManager::getInstance();

  inputPreparation(method_family, program);

  // Generate Calculator
  std::shared_ptr<Scine::Core::Calculator> calc;
  try {
    calc = manager.get<Scine::Core::Calculator>(Scine::Core::Calculator::supports(method_family), program);
  }
  catch (...) {
    return false;
  }
  return true;
}

std::shared_ptr<Scine::Core::Calculator> getCalculator(std::string method_family, std::string program) {
  // Load module manager
  auto& manager = Scine::Core::ModuleManager::getInstance();

  inputPreparation(method_family, program);

  auto split = [&](const std::string& s, char delimiter) -> std::vector<std::string> {
    std::vector<std::string> elements;
    std::stringstream ss;
    ss.str(s);
    std::string item;

    while (std::getline(ss, item, delimiter)) {
      elements.push_back(item);
    }
    return elements;
  };

  // Generate calculator
  const char sep = '/';
  auto listOfMethods = split(method_family, sep);
  auto listOfPrograms = split(program, sep);

  auto dispatcher = [&listOfMethods, &listOfPrograms, &manager]() -> std::vector<std::shared_ptr<Scine::Core::Calculator>> {
    if (listOfMethods.size() != listOfPrograms.size()) {
      throw std::runtime_error("Unequal number of method families and programs given. Please provide a corresponding "
                               "program for each method family in the embedding calculation.");
    }
    if (listOfMethods.size() >= 2 && listOfPrograms.size() >= 2) {
      std::vector<std::shared_ptr<Scine::Core::Calculator>> underlyingCalculators;
      for (long unsigned int i = 0; i < listOfPrograms.size(); i++) {
        underlyingCalculators.emplace_back(manager.get<Scine::Core::Calculator>(
            Scine::Core::Calculator::supports(listOfMethods.at(i)), listOfPrograms.at(i)));
      }
      return underlyingCalculators;
    }
    else {
      throw std::runtime_error(
          "Please provide at least two method families (and the corresponding programs) for an embedding calculation.");
    }
  };

  std::shared_ptr<Scine::Core::Calculator> calc;
  if (listOfMethods.size() < 2 && listOfPrograms.size() < 2) {
    std::shared_ptr<Scine::Core::Calculator> calc;
    try {
      calc = manager.get<Scine::Core::Calculator>(Scine::Core::Calculator::supports(method_family), program);
    }
    catch (...) {
      if (program.empty() || program == "Any") {
        std::cout << "No SCINE module providing '" << method_family << "' is currently loaded.\n";
        std::cout << "Please add the module to the SCINE_MODULE_PATH\n";
        std::cout << "or load the corresponding Python module in order for it to be accessible.\n";
        throw std::runtime_error("Failed to load method/program.");
      }

      std::cout << "No SCINE module named '" << program << "' providing '" << method_family << "' is currently loaded.\n";
      std::cout << "Please add the module to the SCINE_MODULE_PATH\n";
      std::cout << "or load the corresponding Python module in order for it to be accessible.\n";
      throw std::runtime_error("Failed to load method/program.");
    }
    // Return Calculator
    return calc;
  }
  else {
    calc = manager.get<Scine::Core::Calculator>(Scine::Core::Calculator::supports("QMMM"), "Swoose");
    auto castedCalc = std::dynamic_pointer_cast<Scine::Core::EmbeddingCalculator>(calc);
    if (!castedCalc) {
      throw std::runtime_error("Please specify an embedding calculator.");
    }
    auto listOfCalculators = dispatcher();
    castedCalc->setUnderlyingCalculators(listOfCalculators);
    auto settings = castedCalc->settings().getDescriptorCollection();
    return castedCalc;
  }
}

std::shared_ptr<Scine::Core::Calculator> loadSystem(const std::string& path, std::string method_family,
                                                    const pybind11::kwargs& kwargs) {
  // Convert to ValueCollection and read program
  Scine::Utils::UniversalSettings::ValueCollection collection;
  update(collection, pybind11::dict(kwargs), true);
  std::string program = collection.extract("program", std::string{"any"});
  // Generate Calculator
  auto calc = getCalculator(std::move(method_family), program);
  // Apply settings to Calculator
  calc->settings().merge(collection);
  if (!calc->settings().valid()) {
    calc->settings().throwIncorrectSettings();
  }
  // Set initial structure
  auto readResults = Scine::Utils::ChemicalFileHandler::read(path);
  calc->setStructure(readResults.first);
  return calc;
}

Scine::Utils::Settings getAvailableSettings(std::string method_family, std::string program) {
  auto calc = getCalculator(std::move(method_family), std::move(program));
  return calc->settings();
}

Scine::Utils::PropertyList getPossiblePropertiesByStrings(std::string method_family, std::string program) {
  // Get the correct Calculator
  auto calc = getCalculator(std::move(method_family), std::move(program));
  return calc->possibleProperties();
}

void setRequiredProperties(Scine::Core::Calculator& calculator, const std::vector<Scine::Utils::Property>& list) {
  Scine::Utils::PropertyList properties;
  for (const auto& p : list) {
    properties.addProperty(p);
  }

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

  calculator.def_property(
      "log", [](Scine::Core::Calculator& self) -> Scine::Core::Log& { return self.getLog(); },
      [](Scine::Core::Calculator& self, Scine::Core::Log log) -> void { self.setLog(std::move(log)); },
      pybind11::return_value_policy::reference, "Logger of the calculator.");

  calculator.def("set_required_properties", &Scine::Core::Calculator::setRequiredProperties, "Set the required properties");
  calculator.def("set_required_properties", &setRequiredProperties, "Set the required properties");

  calculator.def("get_possible_properties", &Scine::Core::Calculator::possibleProperties,
                 "Yields a list of properties that the calculator can calculate");

  calculator.def("calculate", &Scine::Core::Calculator::calculate, pybind11::arg("dummy") = std::string{},
                 "Execute the calculation");

  calculator.def("get_results", pybind11::overload_cast<>(&Scine::Core::Calculator::results), "Get any stored results.");
  calculator.def("has_results", &hasResults, "Check if results are present.");

  calculator.def("name", &Scine::Core::Calculator::name, "Yields the name of the calculator");

  // Some static helper functions

  m.def("has_calculator", &hasCalculator, pybind11::arg("method_family"), pybind11::arg("program") = "Any",
        "Checks if a calculator with the given method and the given program is available.");
  m.def("get_calculator", &getCalculator, pybind11::arg("method_family"), pybind11::arg("program") = "Any",
        "Generates a calculator with the given method and from the given program.");
  m.def("load_system", &loadSystem, pybind11::arg("path"), pybind11::arg("method_family"),
        "Loads a single system (xyz-file) into a Calculator with the given method and optional settings. (Deprecated)");
  m.def("load_system_into_calculator", &loadSystem, pybind11::arg("path"), pybind11::arg("method_family"),
        "Loads a single system (xyz-file) into a Calculator with the given method and optional settings.");
  m.def("get_available_settings", &getAvailableSettings, pybind11::arg("method_family"), pybind11::arg("program") = "Any",
        "Gives the available default settings of a Calculator with the given method and from the given program");
  m.def("get_possible_properties", &getPossiblePropertiesByStrings, pybind11::arg("method_family"),
        pybind11::arg("program") = "Any",
        "Lists the available properties of a Calculator with the given method and from the given program.");

  // TODO missing: statesHandler (these need to be bound for them to be useful)
}
