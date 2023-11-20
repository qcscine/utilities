/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Settings.h"
#include <Core/Interfaces/Calculator.h>
#include <Core/Interfaces/CalculatorWithReference.h>
#include <Core/Interfaces/WavefunctionOutputGenerator.h>
#include <Utils/Pybind.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <sstream>

namespace PythonHelper {
using VariantType =
    boost::variant<std::shared_ptr<Scine::Core::Calculator>, std::shared_ptr<Scine::Core::CalculatorWithReference>>;
/**
 * @brief Visitor for the variant type. Just casts to the correct type.
 */
struct CasterVisitor : public boost::static_visitor<std::shared_ptr<Scine::Core::WavefunctionOutputGenerator>> {
  using TypeCastTo = std::shared_ptr<Scine::Core::WavefunctionOutputGenerator>;
  template<typename T>
  auto operator()(const T& toCast) const -> TypeCastTo {
    return std::dynamic_pointer_cast<Scine::Core::WavefunctionOutputGenerator>(toCast);
  }
};
} // namespace PythonHelper
void init_wavefunction_output_generator(pybind11::module& m) {
  pybind11::class_<Scine::Core::WavefunctionOutputGenerator, std::shared_ptr<Scine::Core::WavefunctionOutputGenerator>> wavefunctionOutputGenerator(
      m, "WavefunctionOutputGenerator");

  wavefunctionOutputGenerator.doc() =
      "The WavefunctionOutputGenerator is the abstract base for classes obtaining information "
      "based on a wavefunction, f.i. a molden input file from a converged SCF.";

  /**
   * @brief We need to expose a casting utility as in the python side we cannot horizontally cast.
   */
  m.def("to_wf_generator",
        [](const PythonHelper::VariantType& object) -> std::shared_ptr<Scine::Core::WavefunctionOutputGenerator> {
          return boost::apply_visitor(PythonHelper::CasterVisitor(), object);
        });
  /* NOTE: This is absolutely necessary since Core doesn't STORE the constexpr
   * string anywhere, but this way there should be no missing symbols
   */
  constexpr const char* wfOutGenInterfaceStr = Scine::Core::WavefunctionOutputGenerator::interface;
  wavefunctionOutputGenerator.attr("INTERFACE") = wfOutGenInterfaceStr;

  wavefunctionOutputGenerator.def_property(
      "settings",
      [](Scine::Core::WavefunctionOutputGenerator& self) -> Scine::Utils::Settings& { return self.settings(); },
      [](Scine::Core::WavefunctionOutputGenerator& self, Scine::Utils::Settings settings) {
        self.settings() = std::move(settings);
      },
      pybind11::return_value_policy::reference, "Settings of the wavefunction output generator");

  wavefunctionOutputGenerator.def(
      "wavefunction2file",
      pybind11::overload_cast<const std::string&>(&Scine::Core::WavefunctionOutputGenerator::generateWavefunctionInformation),
      "Outputs the wavefunction. For instance, the molecular orbitals for visualization, in a file.");
  wavefunctionOutputGenerator.def(
      "output_wavefunction",
      [&](Scine::Core::WavefunctionOutputGenerator& self) -> std::string {
        std::stringstream sstream;
        self.generateWavefunctionInformation(sstream);
        return sstream.str();
      },
      "Outputs the wavefunction. For instance, the molecular orbitals for visualization.");
}
