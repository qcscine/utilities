/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Core/Interfaces/Calculator.h>
#include <Core/Interfaces/EmbeddingCalculator.h>
#include <Core/ModuleManager.h>
#include <Utils/CalculatorBasics/CalculationRoutines.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/Geometry.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/Settings.h>
#include <Utils/Technical/CloneInterface.h>
#include <Utils/UniversalSettings/ValueCollection.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <algorithm>

class PyStateHandableObject : public Scine::Core::StateHandableObject {
 public:
  /* Inherit the constructors */
  using Scine::Core::StateHandableObject::StateHandableObject;

  void loadState(std::shared_ptr<Scine::Core::State> state) override {
    PYBIND11_OVERRIDE_PURE(void, StateHandableObject, loadState, state);
  }
  std::shared_ptr<Scine::Core::State> getState() const override {
    PYBIND11_OVERRIDE_PURE(std::shared_ptr<Scine::Core::State>, StateHandableObject, getState, );
  }
};

class PyObjectWithStructure : public Scine::Core::ObjectWithStructure {
 public:
  /* Inherit the constructors */
  using Scine::Core::ObjectWithStructure::ObjectWithStructure;
  void setStructure(const Scine::Utils::AtomCollection& structure) override {
    PYBIND11_OVERRIDE_PURE(void, ObjectWithStructure, setStructure, structure);
  }
  std::unique_ptr<Scine::Utils::AtomCollection> getStructure() const override {
    pybind11::gil_scoped_acquire gil; // Acquire the GIL while in this scope.
    pybind11::function override = pybind11::get_override(this, "getStructure");
    if (override) {          // method is found
      auto obj = override(); // Call the Python function.
      Scine::Utils::AtomCollection collection = obj.cast<Scine::Utils::AtomCollection>();
      return std::make_unique<Scine::Utils::AtomCollection>(collection);
    }
    throw std::runtime_error("Missing overload of 'getStructure' in Python Calculator derivative.");
  }
  void modifyPositions(Scine::Utils::PositionCollection newPositions) override {
    PYBIND11_OVERRIDE_PURE(void, ObjectWithStructure, modifyPositions, newPositions);
  }
  const Scine::Utils::PositionCollection& getPositions() const override {
    PYBIND11_OVERRIDE_PURE(const Scine::Utils::PositionCollection&, ObjectWithStructure, getPositions, );
  }
};

template<class Derived, class... Bases>
class CloneHelper : public Bases... {
 public:
  ~CloneHelper() override = default;
  std::shared_ptr<Derived> clone() const {
    return clone_impl();
  }
  virtual std::shared_ptr<Derived> clone_impl() const = 0;

 private:
  Scine::Core::Calculator* cloneImpl() const override {
    throw std::runtime_error("Should never be called");
  }
};

template<class Derived, class... Bases>
class CloneHelper<Scine::Utils::Abstract<Derived>, Bases...> : public Bases... {
 public:
  ~CloneHelper() override = default;

  std::shared_ptr<Derived> clone() const {
    return clone_impl();
    ;
  }
  std::shared_ptr<Scine::Core::Calculator> clone_impl() const = 0;

 private:
  CloneHelper* cloneImpl() const override = 0;
};

class PyCalculator : public Scine::Core::Calculator {
 public:
  /* Inherit the constructors */
  using Scine::Core::Calculator::Calculator;
  PyCalculator() = default;

  PyCalculator(const Scine::Core::Calculator& c) : Scine::Core::Calculator(c){};

  void setRequiredProperties(const Scine::Utils::PropertyList& requiredProperties) override {
    PYBIND11_OVERRIDE_PURE_NAME(void, Calculator, "_set_required_properties_impl", setRequiredProperties, requiredProperties);
  }
  Scine::Utils::PropertyList getRequiredProperties() const override {
    PYBIND11_OVERRIDE_PURE_NAME(Scine::Utils::PropertyList, Scine::Core::Calculator, "_get_required_properties_impl",
                                getRequiredProperties, );
  }
  Scine::Utils::PropertyList possibleProperties() const override {
    PYBIND11_OVERRIDE_PURE_NAME(Scine::Utils::PropertyList, Scine::Core::Calculator, "_possible_properties_impl",
                                possibleProperties, );
  }
  const Scine::Utils::Results& calculate(std::string description) override {
    PYBIND11_OVERRIDE_PURE_NAME(const Scine::Utils::Results&, Scine::Core::Calculator, "_calculate_impl", calculate, description);
  }
  std::string name() const override {
    PYBIND11_OVERRIDE_PURE_NAME(std::string, Scine::Core::Calculator, "_name_impl", name, );
  }
  Scine::Utils::Settings& settings() override {
    PYBIND11_OVERRIDE_PURE_NAME(Scine::Utils::Settings&, Scine::Core::Calculator, "_settings_impl", settings, );
  }
  const Scine::Utils::Settings& settings() const override {
    PYBIND11_OVERRIDE_PURE_NAME(const Scine::Utils::Settings&, Scine::Core::Calculator, "_settings_impl", settings, );
  }
  Scine::Utils::Results& results() override {
    PYBIND11_OVERRIDE_PURE_NAME(Scine::Utils::Results&, Scine::Core::Calculator, "_results_impl", results, );
  }
  const Scine::Utils::Results& results() const override {
    PYBIND11_OVERRIDE_PURE_NAME(const Scine::Utils::Results&, Scine::Core::Calculator, "_results_impl", results, );
  }
  bool supportsMethodFamily(const std::string& methodFamily) const override {
    PYBIND11_OVERRIDE_PURE_NAME(bool, Scine::Core::Calculator, "_supports_method_family_impl", supportsMethodFamily,
                                methodFamily);
  }
  void setStructure(const Scine::Utils::AtomCollection& structure) override {
    PYBIND11_OVERRIDE_PURE_NAME(void, Scine::Core::Calculator, "_set_structure_impl", setStructure, structure);
  }
  std::unique_ptr<Scine::Utils::AtomCollection> getStructure() const override {
    pybind11::gil_scoped_acquire gil; // Acquire the GIL while in this scope.
    // pybind11::function override = pybind11::get_override(this, "getStructure");
    pybind11::function override = pybind11::get_override(this, "_get_structure_impl");
    if (override) {          // method is found
      auto obj = override(); // Call the Python function.
      Scine::Utils::AtomCollection collection = obj.cast<Scine::Utils::AtomCollection>();
      return std::make_unique<Scine::Utils::AtomCollection>(collection);
    }
    throw std::runtime_error("Missing overload of '_get_structure_impl' in Python Calculator derivative.");
  }
  void modifyPositions(Scine::Utils::PositionCollection newPositions) override {
    PYBIND11_OVERRIDE_PURE_NAME(void, Scine::Core::Calculator, "_modify_positions_impl", modifyPositions, newPositions);
  }
  const Scine::Utils::PositionCollection& getPositions() const override {
    PYBIND11_OVERRIDE_PURE_NAME(const Scine::Utils::PositionCollection&, Scine::Core::Calculator, "_get_positions_impl",
                                getPositions, );
  }
  void loadState(std::shared_ptr<Scine::Core::State> state) override {
    PYBIND11_OVERRIDE_PURE_NAME(void, Scine::Core::Calculator, "_load_state_impl", loadState, state);
  }
  std::shared_ptr<Scine::Core::State> getState() const override {
    PYBIND11_OVERRIDE_PURE_NAME(std::shared_ptr<Scine::Core::State>, Scine::Core::Calculator, "_get_state_impl", getState, );
  }
  bool allowsPythonGILRelease() const override {
    PYBIND11_OVERRIDE_PURE_NAME(bool, Scine::Core::Calculator, "_allows_python_gil_release_impl", allowsPythonGILRelease, );
  }

  pybind11::object _clone_impl() const {
    PYBIND11_OVERRIDE_PURE(pybind11::object, Scine::Core::Calculator, _clone_impl, );
  }

  std::shared_ptr<Scine::Core::Calculator> cloneImpl() const override {
    pybind11::gil_scoped_acquire gil; // Acquire the GIL while in this scope.
    pybind11::function override = pybind11::get_override(this, "_clone_impl");
    if (override) { // method is found
      auto cloned = override();
      auto keep_python_state_alive = std::make_shared<pybind11::object>(cloned);
      auto ptr = cloned.cast<Scine::Core::Calculator*>();
      // aliasing shared_ptr: points to `Scine::Core::Calculator* ptr` but refcounts the Python object
      return std::shared_ptr<Scine::Core::Calculator>(keep_python_state_alive, ptr);
    }
    else {
      throw std::runtime_error("Missing overload of '_clone_impl' in Python Calculator derivative.");
    }
  }
};

// Defined in ValueCollectionPython.cpp
extern void update(Scine::Utils::UniversalSettings::ValueCollection& coll, const pybind11::dict& dictionary, bool preserveTypes);

bool hasCalculator(std::string method_family, std::string program) {
  // Load module manager
  auto& manager = Scine::Core::ModuleManager::getInstance();

  Scine::Utils::CalculationRoutines::inputPreparation(method_family, program);

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

std::shared_ptr<Scine::Core::Calculator> loadSystem(const std::string& path, std::string method_family,
                                                    const pybind11::kwargs& kwargs) {
  // Convert to ValueCollection and read program
  Scine::Utils::UniversalSettings::ValueCollection collection;
  update(collection, pybind11::dict(kwargs), true);
  std::string program = collection.extract("program", std::string{"any"});
  // Generate Calculator
  auto calc = Scine::Utils::CalculationRoutines::getCalculator(std::move(method_family), program);
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
  auto calc = Scine::Utils::CalculationRoutines::getCalculator(std::move(method_family), std::move(program));
  return calc->settings();
}

Scine::Utils::PropertyList getPossiblePropertiesByStrings(std::string method_family, std::string program) {
  // Get the correct Calculator
  auto calc = Scine::Utils::CalculationRoutines::getCalculator(std::move(method_family), std::move(program));
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

std::shared_ptr<Scine::Core::Calculator> cloneHelperFunction(const Scine::Core::Calculator& calc) {
  return calc.clone();
}

std::shared_ptr<Scine::Core::Calculator> shared_from_this(Scine::Core::Calculator& calc) {
  return calc.shared_from_this();
}

std::shared_ptr<Scine::Core::EmbeddingCalculator> toEmbedded(Scine::Core::Calculator& calc) {
  auto castedCalc = std::dynamic_pointer_cast<Scine::Core::EmbeddingCalculator>(calc.shared_from_this());
  if (!castedCalc) {
    throw std::runtime_error("Calculator is not an EmbeddingCalculator");
  }
  return castedCalc;
}

void init_calculator(pybind11::module& m) {
  pybind11::class_<Scine::Core::State, std::shared_ptr<Scine::Core::State>> state(m, "State");
  pybind11::class_<Scine::Core::StateHandableObject, PyStateHandableObject, std::shared_ptr<Scine::Core::StateHandableObject>> stateHandableObject(
      m, "StateHandableObject");
  pybind11::class_<Scine::Core::ObjectWithStructure, PyObjectWithStructure, std::shared_ptr<Scine::Core::ObjectWithStructure>> objectWithStructure(
      m, "ObjectWithStructure");
  pybind11::class_<Scine::Core::Calculator, Scine::Core::StateHandableObject, Scine::Core::ObjectWithStructure,
                   PyCalculator, std::shared_ptr<Scine::Core::Calculator>>
      calculator(m, "Calculator");
  pybind11::class_<Scine::Core::EmbeddingCalculator, Scine::Core::Calculator, std::shared_ptr<Scine::Core::EmbeddingCalculator>> embedding_calculator(
      m, "EmbeddingCalculator");
  calculator.def(pybind11::init_alias<>(), "Default Constructor");

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

  calculator.def("get_required_properties", &Scine::Core::Calculator::getRequiredProperties, "Get the required properties");
  calculator.def("set_required_properties", &Scine::Core::Calculator::setRequiredProperties, "Set the required properties");
  calculator.def("set_required_properties", &setRequiredProperties, "Set the required properties");

  calculator.def("get_possible_properties", &Scine::Core::Calculator::possibleProperties,
                 "Yields a list of properties that the calculator can calculate");

  calculator.def("calculate", &Scine::Core::Calculator::calculate, pybind11::arg("dummy") = std::string{},
                 "Execute the calculation");

  calculator.def("get_results", pybind11::overload_cast<>(&Scine::Core::Calculator::results), "Get any stored results.");
  calculator.def("has_results", &hasResults, "Check if results are present.");
  calculator.def(
      "delete_results", [](Scine::Core::Calculator& self) -> void { self.results() = Scine::Utils::Results(); },
      "Deletes the present results.");

  calculator.def("name", &Scine::Core::Calculator::name, "Yields the name of the calculator");
  calculator.def("clone", &cloneHelperFunction, "Yields a copy of the calculator");
  calculator.def("shared_from_this", &shared_from_this, "Yields a shared pointer to copy of the calculator");

  calculator.def("to_embedded_calculator", &toEmbedded,
                 "Tries to transform to an EmbeddingCalculator, raises exception otherwise");
  embedding_calculator.def("get_underlying_calculators", &Scine::Core::EmbeddingCalculator::getUnderlyingCalculators,
                           "Get the underlying calculators");

  // Some static helper functions

  m.def("has_calculator", &hasCalculator, pybind11::arg("method_family"), pybind11::arg("program") = "Any",
        "Checks if a calculator with the given method and the given program is available.");
  m.def("get_calculator", &Scine::Utils::CalculationRoutines::getCalculator, pybind11::arg("method_family"),
        pybind11::arg("program") = "Any", "Generates a calculator with the given method and from the given program.");
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
