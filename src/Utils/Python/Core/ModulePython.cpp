/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Core/Exceptions.h>
#include <Core/Interfaces/Calculator.h>
#include <Core/Module.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>
#include <vector>

class PyModule : public Scine::Core::Module {
 public:
  /* Inherit the constructors */
  using Scine::Core::Module::Module;

  std::string name() const noexcept override {
    PYBIND11_OVERRIDE_PURE(std::string, Module, name, );
  }

  boost::any get(const std::string& interface, const std::string& model) const override {
    pybind11::gil_scoped_acquire gil; // Acquire the GIL while in this scope.
    pybind11::function override = pybind11::get_override(this, "get");
    if (override) { // method is found
      auto obj = override(interface, model);
      if (pybind11::isinstance<Scine::Core::Calculator>(obj)) {
        auto keep_python_state_alive = std::make_shared<pybind11::object>(obj);
        auto ptr = obj.cast<Scine::Core::Calculator*>();
        // aliasing shared_ptr: points to `Scine::Core::Calculator* ptr` but refcounts the Python object
        return std::shared_ptr<Scine::Core::Calculator>(keep_python_state_alive, ptr);
      }
      throw Scine::Core::ClassNotImplementedError();
    }
    throw Scine::Core::ClassNotImplementedError();
  }

  bool has(const std::string& interface, const std::string& model) const noexcept override {
    PYBIND11_OVERRIDE_PURE(bool, Module, has, interface, model);
  }
  std::vector<std::string> announceInterfaces() const noexcept override {
    PYBIND11_OVERRIDE_PURE(std::vector<std::string>, Module, announceInterfaces, );
  }
  std::vector<std::string> announceModels(const std::string& interface) const noexcept override {
    PYBIND11_OVERRIDE_PURE(std::vector<std::string>, Module, announceModels, interface);
  }
};

void init_module(pybind11::module& m) {
  pybind11::class_<Scine::Core::Module, std::shared_ptr<Scine::Core::Module>, PyModule> module(m, "Module");
  module.def(pybind11::init<>());
  module.doc() =
      "Abstract base class for a module, which flexibly provides consumers with common interface derived classes.";
}
