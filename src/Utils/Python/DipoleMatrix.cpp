/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/DataStructures/DipoleMatrix.h"
#include "Utils/Pybind.h"

void init_dipole_matrix(pybind11::module& m) {
  using namespace Scine::Utils;
  pybind11::class_<DipoleMatrix> dipole_matrix(m, "DipoleMatrix");

  /* NOTE: This is a stub just so the matrix is defined for Results and because
   * binding MatrixWithDerivatives (or any methods associated with it) would be
   * a pain.
   */

  dipole_matrix.def("__len__", [](const DipoleMatrix & /* mat */) -> unsigned { return 3; });

  dipole_matrix.def("__getitem__",
                    [](const DipoleMatrix& mat, const unsigned i) -> const Eigen::MatrixXd& { return mat[i]; });
}
