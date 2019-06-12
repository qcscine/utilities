/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Core/Interfaces/Calculator.h>
#include <Core/ModuleManager.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <boost/optional.hpp>
#include <boost/variant.hpp>

using namespace Scine::Core;

// Type caster for boost variant
namespace pybind11 {
namespace detail {

template<typename T>
struct type_caster<boost::optional<T>> : optional_caster<boost::optional<T>> {};

template<typename... Ts>
struct type_caster<boost::variant<Ts...>> : variant_caster<boost::variant<Ts...>> {};

// Specifies the function used to visit the variant -- `apply_visitor` instead of `visit`
template<>
struct visit_helper<boost::variant> {
  template<typename... Args>
  static auto call(Args&&... args) -> decltype(boost::apply_visitor(args...)) {
    return boost::apply_visitor(args...);
  }
};

} // namespace detail
} // namespace pybind11

void module_manager_load(ModuleManager& manager, const std::string& filename) {
  manager.load(filename);
}

using GetVariantType = boost::variant<std::shared_ptr<Calculator>>;

using GetReturnType = boost::optional<GetVariantType>;

GetReturnType module_manager_get(ModuleManager& manager, const std::string& interface, const std::string& model) {
  if (manager.has(interface, model)) {
    if (interface == Calculator::interface) {
      return GetVariantType{manager.get<Calculator>(model)};
    }
  }

  return boost::none;
}

std::vector<std::string> module_manager_models(ModuleManager& manager, const std::string& interface) {
  return manager.getLoadedModels(interface);
}

bool module_manager_has(ModuleManager& manager, const std::string& interface, const std::string& model) {
  return manager.has(interface, model);
}

void init_module_manager(pybind11::module& m) {
  pybind11::class_<ModuleManager, std::unique_ptr<ModuleManager, pybind11::nodelete>> module_manager(m,
                                                                                                     "ModuleManager");
  module_manager.def(
      pybind11::init([]() { return std::unique_ptr<ModuleManager, pybind11::nodelete>(&ModuleManager::getInstance()); }));

  module_manager.def("load", &module_manager_load, "Load a module file");

  module_manager.def_property_readonly("modules", &ModuleManager::getLoadedModuleNames, "List of loaded module names");

  module_manager.def_property_readonly("interfaces", &ModuleManager::getLoadedInterfaces,
                                       "List of interfaces for which at least one model is loaded");

  module_manager.def("models", &module_manager_models, "List of available models of an interface");

  module_manager.def("has", &module_manager_has, "Check whether a particular model of an interface is available");

  module_manager.def("module_loaded", &ModuleManager::moduleLoaded, "Check whether a particular module is loaded");

  module_manager.def("get", &module_manager_get, "Get an instance of an interface model");
}
