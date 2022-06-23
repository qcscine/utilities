/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Math/MachineLearning/Regression/GaussianProcessRegression.h>
#include <Utils/Math/MachineLearning/Regression/KernelRidgeRegression.h>
#include <Utils/Pybind.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils::MachineLearning;

void init_regression_functionalities(pybind11::module& m) {
  auto ml = m.def_submodule("ml");
  ml.doc() = R"delim(
    "Machine Learning Regression Model."
  )delim";

  pybind11::class_<GaussianProcessRegression> gpr{ml, "GaussianProcessRegression", R"delim(
      A class representing a Gaussian process regression.
    )delim"};

  pybind11::class_<GaussianProcessRegression::HyperparameterSpecifier>(ml, "HyperparameterSpecifier")
      .def(pybind11::init<>())
      .def_readwrite("guess", &GaussianProcessRegression::HyperparameterSpecifier::guess)
      .def_readwrite("toOptimize", &GaussianProcessRegression::HyperparameterSpecifier::toOptimize)
      .def_readwrite("bounds", &GaussianProcessRegression::HyperparameterSpecifier::bounds);

  gpr.def(pybind11::init<>());

  gpr.def("train_model", &GaussianProcessRegression::trainModel, pybind11::arg("featureValues"),
          pybind11::arg("targetValues"), R"delim(Trains the GPR model.)delim");

  gpr.def("predict", &GaussianProcessRegression::predict, pybind11::arg("data"), R"delim(Predicts with the GPR model.)delim");

  gpr.def("get_optimized_hyperparams", &GaussianProcessRegression::getOptimizedHyperparameters,
          R"delim(Predicts with the GPR model.)delim");

  gpr.def("get_variance_of_prediction", &GaussianProcessRegression::getVarianceOfPrediction,
          R"delim(Predicts with the GPR model.)delim");

  gpr.def("set_sigmaYSq_guess", &GaussianProcessRegression::setSigmaYSqHyperparametersGuess);

  gpr.def("set_sigmaFSq_guess", &GaussianProcessRegression::setSigmaFSqHyperparametersGuess);

  gpr.def("set_theta_guess", &GaussianProcessRegression::setThetaHyperparametersGuess);

  pybind11::class_<KernelRidgeRegression> krr{ml, "KernelRidgeRegression", R"delim(
      A class representing a kernel ridge regression.
    )delim"};

  krr.def(pybind11::init<>());

  krr.def("train_model", &KernelRidgeRegression::trainModel, pybind11::arg("featureValues"),
          pybind11::arg("targetValues"), R"delim(Trains the kernel ridge regression model.)delim");

  krr.def("predict", &KernelRidgeRegression::predict, pybind11::arg("data"),
          R"delim(Predicts with the kernel ridge regression model.)delim");
}