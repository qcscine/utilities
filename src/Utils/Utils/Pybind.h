/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef INCLUDE_UTILS_PYBIND_H
#define INCLUDE_UTILS_PYBIND_H

#include <Utils/Typenames.h>
#include <pybind11/eigen.h>
#include <pybind11/eval.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <boost/optional.hpp>
#include <boost/variant.hpp>

namespace pybind11 {
namespace detail {

template<typename T>
struct type_caster<boost::optional<T>> : optional_caster<boost::optional<T>> {};

template<typename... Ts>
struct type_caster<boost::variant<Ts...>> : variant_caster<boost::variant<Ts...>> {};

template<>
struct visit_helper<boost::variant> {
  template<typename... Args>
  static auto call(Args&&... args) -> decltype(boost::apply_visitor(args...)) {
    return boost::apply_visitor(args...);
  }
};

} // namespace detail
} // namespace pybind11

namespace Scine {
namespace Utils {

template<typename ContainerHandle, typename StringifierUnary>
inline std::string join(ContainerHandle container, StringifierUnary&& f, std::string sep = ", ") {
  std::string joined;
  for (const auto& element : container) {
    joined += f(element) + sep;
  }
  if (!joined.empty()) {
    joined.erase(joined.size() - sep.size());
  }
  return joined;
}

inline std::string qualifiedName(pybind11::handle handle) {
  auto typehandle = pybind11::type::handle_of(handle);
  auto qualifiedName =
      (typehandle.attr("__module__").cast<std::string>() + "." + typehandle.attr("__qualname__").cast<std::string>());
  return qualifiedName;
}

inline bool isBoundType(pybind11::handle handle) {
  return pybind11::detail::get_internals().registered_types_py.count(Py_TYPE(handle.ptr())) > 0;
}

/**
 * @brief Argument constructor with syntactically valid default argument docstrings
 */
struct Arg {
  Arg(std::string name_) : name(std::move(name_)) {
  }

  // Generate a literal expression for default function argument values
  static std::string compose_literal(pybind11::handle h) {
    namespace py = pybind11;
    auto typehandle = py::type::handle_of(h);
    if (isBoundType(h)) {
      if (py::hasattr(typehandle, "__members__") && py::hasattr(h, "name")) {
        // Bound enum type, can be fully represented
        auto descr = typehandle.attr("__module__").cast<std::string>();
        descr += "." + typehandle.attr("__qualname__").cast<std::string>();
        auto name = h.attr("name");
        if (py::isinstance<py::str>(name)) {
          descr += "." + name.cast<std::string>();
        }
        else {
          descr += "." + name().cast<std::string>();
        }
        return descr;
      }

      // Use ellipsis expression instead of repr to ensure syntactic validity
      return "...";
    }

    if (py::isinstance<py::dict>(h)) {
      std::string literal = "{";
      literal += join(py::reinterpret_borrow<py::dict>(h),
                      [](auto&& v) { return compose_literal(v.first) + ": " + compose_literal(v.second); });
      literal += "}";
      return literal;
    }

    if (py::isinstance<py::list>(h)) {
      std::string literal = "[";
      literal += join(py::reinterpret_borrow<py::list>(h), &compose_literal);
      literal += "]";
      return literal;
    }

    if (py::isinstance<py::tuple>(h)) {
      std::string literal = "(";
      literal += join(py::reinterpret_borrow<py::tuple>(h), &compose_literal);
      literal += ")";
      return literal;
    }

    if (py::isinstance<py::set>(h)) {
      auto v = py::reinterpret_borrow<py::set>(h);
      if (v.empty()) {
        return "set()";
      }
      std::string literal = "{";
      literal += join(v, &compose_literal);
      literal += "}";
      return literal;
    }

    // All other types should be terminal and well-represented by repr
    return repr(h).cast<std::string>();
  }

  template<typename T>
  static pybind11::object make_object(T&& t) {
    return pybind11::reinterpret_steal<pybind11::object>(
        pybind11::detail::make_caster<T>::cast(t, pybind11::return_value_policy::automatic, {}));
  }

  template<typename T>
  pybind11::arg_v operator=(T&& t) const {
    auto object = make_object(std::forward<T>(t));
    std::string description = compose_literal(object);
    return pybind11::arg_v(strdup(name.c_str()), object, strdup(description.c_str()));
  }

  std::string name;
};

template<>
inline pybind11::arg_v Arg::operator=(Position&& position) const {
  auto object = make_object(std::forward<Position>(position));
  std::string description = "numpy." + compose_literal(object);
  return pybind11::arg_v(strdup(name.c_str()), object, strdup(description.c_str()));
}

inline void raiseDeprecationWarning(const std::string& warning) {
  std::string command = "import warnings\n";
  command += "warnings.warn(\"" + warning + "\", DeprecationWarning)";
  pybind11::exec(command);
}

// free function
template<typename Return, typename... Args>
inline auto deprecated(Return (*f)(Args...), const std::string& warning) {
  return [f, warning](Args&&... args) -> decltype(auto) {
    raiseDeprecationWarning(warning);
    return f(args...);
  };
}

// class method
template<typename Return, typename Class, typename... Args>
inline auto deprecated(Return (Class::*f)(Args...), const std::string& warning) {
  return [f, warning](Class* c, Args&&... args) -> decltype(auto) {
    raiseDeprecationWarning(warning);
    return (c->*f)(args...);
  };
}

// const class method
template<typename Return, typename Class, typename... Args>
inline auto deprecated(Return (Class::*f)(Args...) const, const std::string& warning) {
  return [f, warning](const Class* c, Args&&... args) -> decltype(auto) {
    raiseDeprecationWarning(warning);
    return (c->*f)(args...);
  };
}

} // namespace Utils
} // namespace Scine

#endif
