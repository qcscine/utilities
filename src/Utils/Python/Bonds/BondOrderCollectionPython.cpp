/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Bonds/BondOrderCollection.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;

namespace {

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

} // namespace

void init_bond_order_collection(pybind11::module& m) {
  pybind11::class_<BondOrderCollection> bond_order_collection(m, "BondOrderCollection",
                                                              R"delim(
      A sparse symmetric matrix of floating point bond orders

      >>> bo = BondOrderCollection(4)
      >>> bo.empty()
      True
      >>> bo.get_order(1, 2)
      0.0
      >>> bo.set_order(1, 2, 1.3) # sets both (1, 2) and (2, 1)
      >>> bo.empty()
      False
      >>> bo.get_order(1, 2)
      1.3
      >>> bo.get_order(2, 1)
      1.3
      >>> bo.resize(5) # Resizing is non-conservative
      >>> bo.empty()
      True
    )delim");
  bond_order_collection.def(pybind11::init<>());
  bond_order_collection.def(pybind11::init<int>(), "Initialize a bond order collection to a particular size");
  bond_order_collection.def_property("matrix", &BondOrderCollection::getMatrix, &setBOMatrix, "Underlying matrix");
  bond_order_collection.def("resize", &BondOrderCollection::resize, pybind11::arg("N"), R"delim(
    Resize the matrix

    Resizing the matrix does not preserve any existing values
  )delim");
  bond_order_collection.def("set_zero", &BondOrderCollection::setZero, "Set all bond orders to zero");
  bond_order_collection.def("get_system_size", &getBOSystemSize, "Get the system size");
  bond_order_collection.def("set_order", &setBOOrder, pybind11::arg("i"), pybind11::arg("j"), pybind11::arg("order"),
                            R"delim(
      Set a bond order

      Sets both (i, j) and (j, i) entries in the matrix

      :param i: First matrix index
      :param j: Second matrix index
      :param order: Order to set
    )delim");
  bond_order_collection.def("get_order", &getBOOrder, pybind11::arg("i"), pybind11::arg("j"), "Get a bond order");
  bond_order_collection.def("empty", &BondOrderCollection::empty, "Checks whether there are no entries in the BO matrix");

  // Operators
  bond_order_collection.def(pybind11::self == pybind11::self);
  bond_order_collection.def(pybind11::self != pybind11::self);
}
