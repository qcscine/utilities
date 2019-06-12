/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Bonds/BondOrderCollection.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;

namespace detail {

void setBOMatrix(BondOrderCollection& bo, Eigen::SparseMatrix<double> matr) {
  bo.setMatrix(std::move(matr));
}

int getBOSystemSize(const BondOrderCollection& bo) {
  return bo.getSystemSize();
}

void setBOOrder(BondOrderCollection& bo, int i, int j, double order) {
  bo.setOrder(i, j, order);
}

double getBOOrder(const BondOrderCollection& bo, int i, int j) {
  return bo.getOrder(i, j);
}

} // namespace detail

void init_bond_order_collection(pybind11::module& m) {
  pybind11::class_<BondOrderCollection> bond_order_collection(m, "BondOrderCollection");
  bond_order_collection.def(pybind11::init<int>());
  bond_order_collection.def_property("matrix", &BondOrderCollection::getMatrix, &detail::setBOMatrix);
  bond_order_collection.def("resize", &BondOrderCollection::resize, "Resize the matrix");
  bond_order_collection.def("set_zero", &BondOrderCollection::setZero, "Set all bond orders to zero");
  bond_order_collection.def("get_system_size", &detail::getBOSystemSize, "Get the system size");
  bond_order_collection.def("set_order", &detail::setBOOrder, "Set a bond order");
  bond_order_collection.def("get_order", &detail::getBOOrder, "Get a bond order");
  bond_order_collection.def("empty", &BondOrderCollection::empty, "Checks whether there are no entries in the BO matrix");

  // Operators
  bond_order_collection.def(pybind11::self == pybind11::self);
  bond_order_collection.def(pybind11::self != pybind11::self);
}
