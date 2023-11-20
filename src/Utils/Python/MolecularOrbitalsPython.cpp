/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/DataStructures/MolecularOrbitals.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <boost/optional.hpp>

using namespace Scine::Utils;

void init_molecular_orbitals(pybind11::module& m) {
  pybind11::class_<MolecularOrbitals> molecular_orbitals(m, "MolecularOrbitals");

  molecular_orbitals.def_property_readonly(
      "alpha_matrix", pybind11::overload_cast<>(&MolecularOrbitals::alphaMatrix, pybind11::const_),
      "Returns the coefficient matrix for electrons with alpha spin from an unrestricted calculation.");
  molecular_orbitals.def_property_readonly(
      "beta_matrix", pybind11::overload_cast<>(&MolecularOrbitals::betaMatrix, pybind11::const_),
      "Returns the coefficient matrix for electrons with beta spin from an unrestricted calculation.");
  molecular_orbitals.def_property_readonly(
      "restricted_matrix", pybind11::overload_cast<>(&MolecularOrbitals::restrictedMatrix, pybind11::const_),
      "Returns the coefficient matrix for all electrons from a restricted calculation.");
}
