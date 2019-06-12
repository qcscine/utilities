/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/CalculatorBasics/Results.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <boost/optional.hpp>

// boost optional type caster
namespace pybind11 {
namespace detail {
template<typename T>
struct type_caster<boost::optional<T>> : optional_caster<boost::optional<T>> {};
} // namespace detail
} // namespace pybind11

using namespace Scine::Utils;

/* Adds a property that might be None to the python class */
template<typename T, class HasMethod, class GetMethod, class SetMethod>
void add_optional_property(pybind11::class_<Results>& pyclass_results, const char* method_name, HasMethod&& has_method,
                           GetMethod&& get_method, SetMethod&& set_method) {
  pyclass_results.def_property(method_name,
                               [=](Results& results) -> boost::optional<T> {
                                 if (!std::mem_fn(has_method)(results)) {
                                   return boost::none;
                                 }

                                 return std::mem_fn(get_method)(results);
                               },
                               [=](Results& results, const boost::optional<T>& valueOptional) -> void {
                                 if (valueOptional) {
                                   std::mem_fn(set_method)(results, *valueOptional);
                                 }
                               });
}

void init_results(pybind11::module& m) {
  pybind11::class_<Results> results(m, "Results");

  results.def(pybind11::init<>());

  add_optional_property<std::string>(results, "description", &Results::hasDescription, &Results::getDescription,
                                     &Results::setDescription);

  add_optional_property<double>(results, "energy", &Results::hasEnergy, &Results::getEnergy, &Results::setEnergy);

  add_optional_property<GradientCollection>(results, "gradients", &Results::hasGradients, &Results::getGradients,
                                            &Results::setGradients);

  add_optional_property<HessianMatrix>(results, "hessian", &Results::hasHessian, &Results::getHessian, &Results::setHessian);

  add_optional_property<Dipole>(results, "dipole", &Results::hasDipole, &Results::getDipole, &Results::setDipole);

  add_optional_property<DipoleGradient>(results, "dipole_gradient", &Results::hasDipoleGradient,
                                        &Results::getDipoleGradient, &Results::setDipoleGradient);

  add_optional_property<DipoleMatrix>(results, "ao_dipole_matrix", &Results::hasAODipoleMatrix,
                                      &Results::getAODipoleMatrix, &Results::setAODipoleMatrix);

  add_optional_property<DipoleMatrix>(results, "mo_dipole_matrix", &Results::hasMODipoleMatrix,
                                      &Results::getMODipoleMatrix, &Results::setMODipoleMatrix);

  add_optional_property<Eigen::MatrixXd>(results, "one_electron_matrix", &Results::hasOneElectronMatrix,
                                         &Results::getOneElectronMatrix, &Results::setOneElectronMatrix);

  add_optional_property<SpinAdaptedMatrix>(results, "two_electron_matrix", &Results::hasTwoElectronMatrix,
                                           &Results::getTwoElectronMatrix, &Results::setTwoElectronMatrix);

  add_optional_property<BondOrderCollection>(results, "bond_orders", &Results::hasBondOrders, &Results::getBondOrders,
                                             &Results::setBondOrders);
}
