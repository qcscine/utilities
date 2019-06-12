/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Geometry/StructuralCompletion.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;

void init_structural_completion(pybind11::module& m) {
  pybind11::class_<StructuralCompletion> structural_completion(m, "StructuralCompletion");

  structural_completion.def_static("generate_three_tetrahedron_corners_from_one",
                                   &StructuralCompletion::generate3TetrahedronCornersFrom1Other,
                                   "Generates three missing positions in a tetrahedron");
  structural_completion.def_static("generate_two_tetrahedron_corners_from_two",
                                   &StructuralCompletion::generate2TetrahedronCornersFrom2Others,
                                   "Generates two missing positions in a tetrahedron");
  structural_completion.def_static("generate_one_tetrahedron_corner_from_three",
                                   &StructuralCompletion::generate1TetrahedronCornerFrom3Others,
                                   "Generates one missing position in a tetrahedron");
  structural_completion.def_static("generate_one_triangle_corner_from_two",
                                   &StructuralCompletion::generate1TriangleCornerFrom2Others,
                                   "Generates one triangle corner positions from two");
  structural_completion.def_static("generate_two_triangle_corners_from_one",
                                   &StructuralCompletion::generate2TriangleCornersFrom1Other,
                                   "Generates two triangle corner positions from one");
}
