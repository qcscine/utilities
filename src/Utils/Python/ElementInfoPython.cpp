/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Geometry/ElementInfo.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;

void init_element_info(pybind11::module& m) {
  pybind11::class_<ElementInfo> element_info(m, "ElementInfo");

  element_info.def_static("element_from_symbol", &ElementInfo::elementTypeForSymbol, pybind11::arg("element_str"),
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
  element_info.def_static("symbol", &ElementInfo::symbol, pybind11::arg("element"),
                          "Translate an ElementType into its string representation");
  element_info.def_static("all_implemented_elements", &ElementInfo::allImplementedElements,
                          "Gives a list of all implemented ElementTypes");
  element_info.def_static("mass", &ElementInfo::mass, pybind11::arg("element"),
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
  element_info.def_static("covalent_radius", &ElementInfo::covalentRadius, pybind11::arg("element"),
                          R"delim(
      Returns the covalent radius of an element in atomic units.

      References:

      - Atomic Radii of the Elements in CRC Handbook of Chemistry and Physics,
        100th Edition (Internet Version 2019), John R. Rumble, ed., CRC
        Press/Taylor & Francis, Boca Raton, FL.
      - DOI: 10.1039/b801115j
      - DOI: 10.1002/chem.200800987

      :param element: The element type for which to fetch the covalent radius
      :return: covalent radius in atomic units

      >>> ElementInfo.covalent_radius(ElementType.H)
      0.604712360146505
    )delim");
  element_info.def_static("vdw_radius", &ElementInfo::vdwRadius, pybind11::arg("element"),
                          "van der Waals radius in atomic units");
  element_info.def_static("Z", &ElementInfo::Z, pybind11::arg("element"),
                          R"delim(
      Atomic number of an element

      >>> ElementInfo.Z(ElementType.H)
      1
      >>> ElementInfo.Z(ElementType.C14)
      6
    )delim");
  element_info.def_static("A", &ElementInfo::A, pybind11::arg("element"),
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
  element_info.def_static("abundance", &ElementInfo::abundance, pybind11::arg("element"),
                          R"delim(
      Natural abundance of an isotope

      Raises RuntimeError if atomic mass unspecified. The stored natural
      abundances of particular isotopes may not sum to one, but may all be zero
      for cases in which no natural abundances have been measured.

      >>> ElementInfo.abundance(ElementType.H)
      Traceback (most recent call last):
        ...
      RuntimeError: Unspecified isotope has no abundance
      >>> ElementInfo.abundance(ElementType.H1)
      0.999885
      >>> ElementInfo.abundance(ElementType.D)
      0.000115
    )delim");
  element_info.def_static("base", &ElementInfo::base, pybind11::arg("element"),
                          R"delim(
      Returns the base of an isotope (e.g. Li for Li6)

      >>> assert ElementInfo.base(ElementType.Li6) == ElementType.Li
      >>> assert ElementInfo.base(ElementType.H) == ElementType.H
      >>> assert ElementInfo.base(ElementType.H1) == ElementType.H
    )delim");
  element_info.def_static("element", &ElementInfo::element, pybind11::arg("z"),
                          R"delim(
      Compose an element from an atomic number

      >>> ElementInfo.element(6)
      ElementType.C
      >>> ElementInfo.element(4)
      ElementType.Be
    )delim");
  element_info.def_static("isotope", &ElementInfo::isotope, pybind11::arg("z"), pybind11::arg("a"),
                          R"delim(
      Compose an isotope from an atomic number and an atomic mass number

      :param z: The atomic number (number of protons)
      :param a: The atomic mass number (number of protons and neutrons)
      :rtype: ElementType

      >>> assert ElementInfo.isotope(6, 12) == ElementType.C12
    )delim");
  element_info.def_static("isotopes", &ElementInfo::isotopes, pybind11::arg("element"),
                          R"delim(
      Returns isotopes of an element, in unsorted order

      >>> assert sorted(ElementInfo.isotopes(ElementType.H)) == [ElementType.H1, ElementType.D, ElementType.T]
    )delim");
  element_info.def_static("val_electrons", &ElementInfo::valElectrons, pybind11::arg("element"), "Number of valence electrons");
  element_info.def_static("s_electrons", &ElementInfo::sElectrons, pybind11::arg("element"), "Number of s valence electrons");
  element_info.def_static("p_electrons", &ElementInfo::pElectrons, pybind11::arg("element"), "Number of p valence electrons");
  element_info.def_static("d_electrons", &ElementInfo::dElectrons, pybind11::arg("element"), "Number of d valence electrons");
}
