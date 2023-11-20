/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Core/Log.h>
#include <Utils/Math/IterativeDiagonalizer/DavidsonDiagonalizer.h>
#include <Utils/Math/IterativeDiagonalizer/IndirectPreconditionerEvaluator.h>
#include <Utils/Math/IterativeDiagonalizer/IndirectSigmaVectorEvaluator.h>
#include <Utils/Settings.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;

class PySigmaVectorEvaluator : public SigmaVectorEvaluator {
 public:
  using SigmaVectorEvaluator::SigmaVectorEvaluator;

  auto evaluate(const Eigen::MatrixXd& guessVectors) const -> const Eigen::MatrixXd& override {
    PYBIND11_OVERLOAD_PURE(const Eigen::MatrixXd&, SigmaVectorEvaluator, evaluate, guessVectors);
  }
  void collapsed(int newSubspaceDimension) override {
    PYBIND11_OVERLOAD_PURE(void, SigmaVectorEvaluator, collapsed, newSubspaceDimension);
  }
};

class PyPreconditionerEvaluator : public PreconditionerEvaluator {
 public:
  auto evaluate(const Eigen::VectorXd& toPrecondition, double eigenvalue) const -> Eigen::VectorXd override {
    PYBIND11_OVERLOAD_PURE(Eigen::VectorXd, PreconditionerEvaluator, evaluate, toPrecondition, eigenvalue);
  }
};

void init_iterative_diagonalizer(pybind11::module& m) {
  pybind11::class_<SigmaVectorEvaluator, PySigmaVectorEvaluator, std::shared_ptr<SigmaVectorEvaluator>>(
      m, "SigmaVectorEvaluator")
      .def(pybind11::init<>())
      .def("evaluate", &SigmaVectorEvaluator::evaluate)
      .def("collapsed", &SigmaVectorEvaluator::collapsed);

  pybind11::class_<PreconditionerEvaluator, PyPreconditionerEvaluator, std::shared_ptr<PreconditionerEvaluator>>(
      m, "PreconditionerEvaluator")
      .def(pybind11::init<>())
      .def("evaluate", &PreconditionerEvaluator::evaluate);

  pybind11::class_<IndirectSigmaVectorEvaluator<Eigen::MatrixXd>, SigmaVectorEvaluator,
                   std::shared_ptr<IndirectSigmaVectorEvaluator<Eigen::MatrixXd>>>(m, "IndirectSigmaVectorEvaluator")
      .def(pybind11::init([](pybind11::EigenDRef<Eigen::MatrixXd> m) {
        return std::make_shared<IndirectSigmaVectorEvaluator<Eigen::MatrixXd>>(m);
      }))
      .def("evaluate", &IndirectSigmaVectorEvaluator<Eigen::MatrixXd>::evaluate)
      .def("subspaceCollapsed", &IndirectSigmaVectorEvaluator<Eigen::MatrixXd>::collapsed);

  pybind11::class_<IndirectPreconditionerEvaluator, PreconditionerEvaluator, std::shared_ptr<IndirectPreconditionerEvaluator>>(
      m, "IndirectPreconditionerEvaluator")
      .def(pybind11::init<const Eigen::VectorXd&>())
      .def("evaluate", &IndirectPreconditionerEvaluator::evaluate);

  pybind11::class_<NonOrthogonalDavidson> nonOrthogonalDavidson(m, "NonOrthogonalDavidson");
  nonOrthogonalDavidson.def(pybind11::init<int, int>());
  nonOrthogonalDavidson.def_property(
      "sigma_vector_evaluator", pybind11::overload_cast<>(&NonOrthogonalDavidson::getSigmaVectorEvaluator, pybind11::const_),
      [](NonOrthogonalDavidson& self, std::shared_ptr<SigmaVectorEvaluator> eval) {
        self.setSigmaVectorEvaluator(std::move(eval));
      },
      pybind11::return_value_policy::reference, "The sigma vector evaluator.");
  nonOrthogonalDavidson.def("set_preconditioner", &NonOrthogonalDavidson::setPreconditionerEvaluator, "Sets the preconditioner.");
  nonOrthogonalDavidson.def_property(
      "settings", [](NonOrthogonalDavidson& self) -> Settings& { return self.settings(); },
      [](NonOrthogonalDavidson& self, Settings settings) { self.settings() = std::move(settings); },
      pybind11::return_value_policy::reference, "Settings of the Davidson diagonalizer.");

  nonOrthogonalDavidson.def("apply_settings", &NonOrthogonalDavidson::applySettings, "Applies the settings given.");
  nonOrthogonalDavidson.def("set_guess", [](NonOrthogonalDavidson& self, pybind11::EigenDRef<Eigen::MatrixXd> m) {
    self.setGuess(Eigen::MatrixXd(m));
  });
  nonOrthogonalDavidson.def_property_readonly("eigenpairs", &NonOrthogonalDavidson::getEigenPairs,
                                              "The solution of the diagonalization.");
  nonOrthogonalDavidson.def("solve", &NonOrthogonalDavidson::solve,
                            "Solve the diagonalization with the given sigma vector evaluator and preconditioner.");

  pybind11::class_<OrthogonalDavidson> orthogonalDavidson(m, "OrthogonalDavidson");
  orthogonalDavidson.def(pybind11::init<int, int>());
  orthogonalDavidson.def_property(
      "sigma_vector_evaluator", pybind11::overload_cast<>(&OrthogonalDavidson::getSigmaVectorEvaluator, pybind11::const_),
      [](OrthogonalDavidson& self, std::shared_ptr<SigmaVectorEvaluator> eval) {
        self.setSigmaVectorEvaluator(std::move(eval));
      },
      pybind11::return_value_policy::reference, "The sigma vector evaluator.");
  orthogonalDavidson.def("set_preconditioner", &OrthogonalDavidson::setPreconditionerEvaluator, "Sets the preconditioner.");
  orthogonalDavidson.def_property(
      "settings", [](OrthogonalDavidson& self) -> Settings& { return self.settings(); },
      [](OrthogonalDavidson& self, Settings settings) { self.settings() = std::move(settings); },
      pybind11::return_value_policy::reference, "Settings of the Davidson diagonalizer.");

  orthogonalDavidson.def("apply_settings", &OrthogonalDavidson::applySettings, "Applies the settings given.");
  orthogonalDavidson.def("set_guess", [](OrthogonalDavidson& self, pybind11::EigenDRef<Eigen::MatrixXd> m) {
    self.setGuess(Eigen::MatrixXd(m));
  });
  orthogonalDavidson.def_property_readonly("eigenpairs", &OrthogonalDavidson::getEigenPairs,
                                           "The solution of the diagonalization.");
  orthogonalDavidson.def("solve", &OrthogonalDavidson::solve,
                         "Solve the diagonalization with the given sigma vector evaluator and preconditioner.");
}
