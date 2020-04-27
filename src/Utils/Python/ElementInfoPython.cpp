/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Geometry/ElementInfo.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;

void init_element_info(pybind11::module& m) {
  pybind11::class_<ElementInfo> element_info(m, "ElementInfo");

  element_info.def_static("element_from_symbol", &ElementInfo::elementTypeForSymbol,
                          R"delim(
      Translate a string representation to an ElementType

      Permissive regarding digits specifying isotopic atomic mass numbers, either
      pre- or postfixed.

      >>> assert ElementInfo.element_from_symbol("H") == ElementType.H
      >>> assert ElementInfo.element_from_symbol("H1") == ElementType.H1
      >>> assert ElementInfo.element_from_symbol("1H") == ElementType.H1
      >>> assert ElementInfo.element_from_symbol("D") == ElementType.D
      >>> assert ElementInfo.element_from_symbol("2H") == ElementType.D
      >>> assert ElementInfo.element_from_symbol("T") == ElementType.T
      >>> assert ElementInfo.element_from_symbol("3H") == ElementType.T
    )delim");
  element_info.def_static("symbol", &ElementInfo::symbol, "Translate an ElementType into its string representation");
  element_info.def_static("mass", &ElementInfo::mass,
                          R"delim(
      Standard atomic weight of element type

      The standard atomic weight of an element (e.g. H) is the average of its
      isotopic weights weighted by their natural abundance. If no natural
      abundance for an element was measured or no standard atomic weight is
      defined, returns the weight of one of its isotopes.

      The atomic weight of an isotope (e.g. D) is the mass of the isotope scaled
      onto the standard atomic weight scale, where the standard atomic weight of
      C-12 is set to 12.

      :param element: The element type for which to fetch the standard atomic weight
      :return: standard atomic weight in unified atomic mass units (u)

      >>> ElementInfo.mass(ElementType.H) # H is composed of H1, D and T
      1.0079
      >>> ElementInfo.mass(ElementType.D)
      2.01410177812
    )delim");
  element_info.def_static("vdw_radius", &ElementInfo::vdwRadius, "van der Waals radius in atomic units");
  element_info.def_static("Z", &ElementInfo::Z,
                          R"delim(
      Atomic number of an element

      >>> ElementInfo.Z(ElementType.H)
      1
      >>> ElementInfo.Z(ElementType.C14)
      6
    )delim");
  element_info.def_static("A", &ElementInfo::A,
                          R"delim(
      Atomic mass number of an element

      Returns zero for non-monoisotopic elements.

      >>> ElementInfo.A(ElementType.H) # H has isotopes H1, D and T
      0
      >>> ElementInfo.A(ElementType.D) # For isotopes, A is nonzero
      2
      >>> ElementInfo.A(ElementType.Be) # Be is monoisotopic, A is also nonzero
      9
    )delim");
  element_info.def_static("abundance", &ElementInfo::abundance,
                          R"delim(
      Natural abundance of an isotope

      Returns zero for non-monoisotopic elements. The stored natural abundances
      of particular isotopes may not sum to one, but may all be zero for cases
      in which no natural abundances have been measured.

      >>> ElementInfo.abundance(ElementType.H)
      0.0
      >>> ElementInfo.abundance(ElementType.H1)
      0.999885
    )delim");
  element_info.def_static("base", &ElementInfo::base,
                          R"delim(
      Returns the base of an isotope (e.g. Li for Li6)

      >>> assert ElementInfo.base(ElementType.Li6) == ElementType.Li
      >>> assert ElementInfo.base(ElementType.H) == ElementType.H
      >>> assert ElementInfo.base(ElementType.H1) == ElementType.H
    )delim");
  element_info.def_static("element", &ElementInfo::element,
                          R"delim(
      Compose an element from an atomic number

      >>> ElementInfo.element(6) == ElementType.C
      >>> ElementInfo.element(4) == ElementType.Be
    )delim");
  element_info.def_static("isotope", &ElementInfo::isotope,
                          R"delim(
      Compose an isotope from an atomic number and an atomic mass number

      :param z: The atomic number (number of protons)
      :param a: The atomic mass number (number of protons and neutrons)
      :rtype: ElementType

      >>> assert ElementInfo.isotope(6, 12) == ElementType.C12
    )delim");
  element_info.def_static("isotopes", &ElementInfo::isotopes,
                          R"delim(
      Returns isotopes of an element, in unsorted order

      >>> assert sorted(ElementInfo.isotopes(ElementType.H)) == [ElementType.H1, ElementType.D, ElementType.T]
    )delim");
  element_info.def_static("val_electrons", &ElementInfo::valElectrons, "Number of valence electrons");
  element_info.def_static("s_electrons", &ElementInfo::sElectrons, "Number of s valence electrons");
  element_info.def_static("p_electrons", &ElementInfo::pElectrons, "Number of p valence electrons");
  element_info.def_static("d_electrons", &ElementInfo::dElectrons, "Number of d valence electrons");
}
