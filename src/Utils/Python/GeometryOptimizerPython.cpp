/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/GeometryOptimization/GeometryOptimization.h>
#include <Utils/GeometryOptimization/GeometryOptimizer.h>
#include <Utils/Optimizer/GradientBased/Bfgs.h>
#include <Utils/Optimizer/GradientBased/Dimer.h>
#include <Utils/Optimizer/GradientBased/Lbfgs.h>
#include <Utils/Optimizer/GradientBased/SteepestDescent.h>
#include <Utils/Optimizer/HessianBased/Bofill.h>
#include <Utils/Optimizer/HessianBased/EigenvectorFollowing.h>
#include <Utils/Optimizer/HessianBased/NewtonRaphson.h>
#include <Utils/Pybind.h>
#include <Utils/UniversalSettings/ValueCollection.h>
#include <pybind11/functional.h>

using namespace Scine;
using namespace Utils;

namespace {

enum class OptimizerName { SteepestDescent, Bofill, Dimer, EigenvectorFollowing, Bfgs, Lbfgs, NewtonRaphson };

using GeometryOptimizationObserver = std::function<void(int, double, const AtomCollection& atomCollection)>;

template<typename OptimizerType>
AtomCollection optimizeBase(Core::Calculator& calculator, Core::Log& logger,
                            const GeometryOptimizationObserver& observer, Settings optSettings) {
  GeometryOptimizer<OptimizerType> geometryOptimizer{calculator};

  // unfitting setting for optimizer has to be extracted first
  bool sanityCheck = optSettings.extract(std::string("sanity_check"), true);
  if (optSettings.name() != "empty") {
    Settings settings = geometryOptimizer.getSettings();
    settings.merge(optSettings);
    geometryOptimizer.setSettings(settings);
  }
  if (sanityCheck && !GeometryOptimization::settingsMakeSense(geometryOptimizer)) {
    throw std::logic_error("The given calculator settings are too inaccurate for the given convergence criteria of "
                           "this optimization Task.\n"
                           "If you want to ignore this sanity check, add the setting 'sanity_check': False "
                           "to the settings.");
  }

  if (observer) {
    geometryOptimizer.addObserver([&](const int& cycle, const double& energy, const Eigen::VectorXd & /* parameters */) -> void {
      auto structurePtr = calculator.getStructure();
      observer(cycle, energy, *structurePtr);
    });
  }

  auto structurePtr = calculator.getStructure();
  // Core::Log logger = Core::Log::silent();
  geometryOptimizer.optimize(*structurePtr, logger);
  return *structurePtr;
}

Settings optimization_settings(Core::Calculator& calculator, OptimizerName optimizer = OptimizerName::SteepestDescent) {
  switch (optimizer) {
    case OptimizerName::SteepestDescent:
      return GeometryOptimizer<SteepestDescent>(calculator).getSettings();

    case OptimizerName::Bofill:
      return GeometryOptimizer<Bofill>(calculator).getSettings();

    case OptimizerName::Dimer:
      return GeometryOptimizer<Dimer>(calculator).getSettings();

    case OptimizerName::EigenvectorFollowing:
      return GeometryOptimizer<EigenvectorFollowing>(calculator).getSettings();

    case OptimizerName::Bfgs:
      return GeometryOptimizer<Bfgs>(calculator).getSettings();

    case OptimizerName::Lbfgs:
      return GeometryOptimizer<Lbfgs>(calculator).getSettings();

    case OptimizerName::NewtonRaphson:
      return GeometryOptimizer<NewtonRaphson>(calculator).getSettings();

    default:
      throw std::logic_error("Invalid Optimizer selected.");
  }
}

AtomCollection optimize(Core::Calculator& calculator, Core::Log& logger, OptimizerName optimizer = OptimizerName::SteepestDescent,
                        const GeometryOptimizationObserver& observer = {}, const Settings& settings = Settings("empty")) {
  switch (optimizer) {
    case OptimizerName::SteepestDescent:
      return optimizeBase<SteepestDescent>(calculator, logger, observer, settings);

    case OptimizerName::Bofill:
      return optimizeBase<Bofill>(calculator, logger, observer, settings);

    case OptimizerName::Dimer:
      return optimizeBase<Dimer>(calculator, logger, observer, settings);

    case OptimizerName::EigenvectorFollowing:
      return optimizeBase<EigenvectorFollowing>(calculator, logger, observer, settings);

    case OptimizerName::Bfgs:
      return optimizeBase<Bfgs>(calculator, logger, observer, settings);

    case OptimizerName::Lbfgs:
      return optimizeBase<Lbfgs>(calculator, logger, observer, settings);

    case OptimizerName::NewtonRaphson:
      return optimizeBase<NewtonRaphson>(calculator, logger, observer, settings);

    default:
      throw std::logic_error("Invalid Optimizer selected.");
  }
}

} // namespace
void init_geometry_optimize(pybind11::module& m) {
  pybind11::enum_<OptimizerName> optimizer(m, "Optimizer");
  optimizer.value("SteepestDescent", OptimizerName::SteepestDescent)
      .value("Bofill", OptimizerName::Bofill)
      .value("Dimer", OptimizerName::Dimer)
      .value("EigenvectorFollowing", OptimizerName::EigenvectorFollowing)
      .value("Bfgs", OptimizerName::Bfgs)
      .value("Lbfgs", OptimizerName::Lbfgs)
      .value("NewtonRaphson", OptimizerName::NewtonRaphson);

  m.def("geometry_optimization_settings", &optimization_settings, pybind11::arg("calculator"),
        Arg("optimizer") = OptimizerName::SteepestDescent, "Settings available for geometry optimization");

  m.def("geometry_optimize", &optimize, pybind11::arg("calculator"), pybind11::arg("logger"),
        Arg("optimizer") = OptimizerName::SteepestDescent, pybind11::arg("observer") = GeometryOptimizationObserver{},
        pybind11::arg("settings") = Settings{"empty"},
        R"delim(
    Geometry optimize a structure using a calculator

    :param calculator: Calculator with which to calculate energy, gradient and
      (if required in the optimizer, the hessian). Store the structure you want
      to optimize within the calculator before passing it to this function.
    :param logger: The logger to which eventual output is written.
    :param optimizer: The optimizer with which to perform the geometry
      minimization.
    :param observer: A function of signature (int, float, AtomCollection) -> None
      that is called in each iteration with the current cycle number, energy, and
      structure.
    :param settings: Optional additional settings for the optimization.

    :returns: The optimized structure
  )delim");
}
