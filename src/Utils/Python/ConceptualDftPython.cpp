/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Properties/Reactivity/ConceptualDft.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <boost/optional.hpp>

using namespace Scine::Utils::ConceptualDft;

void init_conceptual_dft(pybind11::module& m) {
  auto conceptual_dft_submodule = m.def_submodule("conceptual_dft", R"delim(
      Functionality to calculate various conceptual DFT properties based on the finite difference approximation.

      To use this functionality the energy and/or atomic charges of the optimized structure of interest and of the same
      structure (not reoptimized!) with one additional electron and one electron less are needed.

      An introduction into conceptual DFT can be found in:
      Chattaraj, P. K. Chemical Reactivity Theory : A Density Functional View; CRC Press, 2009.
      https://doi.org/10.1201/9781420065442.
    )delim");

  pybind11::class_<GlobalConceptualDftContainer> globalConceptualDftContainer(conceptual_dft_submodule,
                                                                              "GlobalConceptualDftContainer");
  globalConceptualDftContainer.def_readonly("chemical_potential", &GlobalConceptualDftContainer::chemicalPotential);
  globalConceptualDftContainer.def_readonly("electronegativity", &GlobalConceptualDftContainer::electronegativity);
  globalConceptualDftContainer.def_readonly("hardness", &GlobalConceptualDftContainer::hardness);
  globalConceptualDftContainer.def_readonly("softness", &GlobalConceptualDftContainer::softness);
  globalConceptualDftContainer.def_readonly("electrophilicity", &GlobalConceptualDftContainer::electrophilicity);

  pybind11::class_<LocalConceptualDftContainer> localConceptualDftContainer(conceptual_dft_submodule,
                                                                            "LocalConceptualDftContainer");
  localConceptualDftContainer.def_readonly("fukui_plus", &LocalConceptualDftContainer::fukuiPlus);
  localConceptualDftContainer.def_readonly("fukui_minus", &LocalConceptualDftContainer::fukuiMinus);
  localConceptualDftContainer.def_readonly("fukui_radical", &LocalConceptualDftContainer::fukuiRadical);
  localConceptualDftContainer.def_readonly("dual_descriptor", &LocalConceptualDftContainer::dualDescriptor);

  pybind11::class_<ConceptualDftContainer> conceptualDftContainer(conceptual_dft_submodule, "ConceptualDftContainer");
  // NOTE: Cannot define global or local, those are reserved keywords in python
  conceptualDftContainer.def_readonly("global_v", &ConceptualDftContainer::global);
  conceptualDftContainer.def_readonly("local_v", &ConceptualDftContainer::local);

  conceptual_dft_submodule.def("calculate", &calculate, pybind11::arg("energy"), pybind11::arg("atomic_charges"),
                               pybind11::arg("energy_plus"), pybind11::arg("atomic_charges_plus"),
                               pybind11::arg("energy_minus"), pybind11::arg("atomic_charges_minus"), R"delim(
      Calculates a set of global and local conceptual DFT parameters.
      Note: The quality of the resulting Fukui and dual descriptor indices heavily depends on the quality of the atomic
      charges.
      We recommend Hirshfeld charges.

      :param energy: The energy of the structure of interest
      :param: atomic_charges The atomic charges of the structure of interest
      :param energy_plus: The energy for the same geometry with one additional electron
      :param: atomic_charges_plus The atomic charges for the same geometry with one additional electron
      :param energy_minus: The energy for the same geometry with one electron less
      :param: atomic_charges_minus The atomic charges for the same geometry with one electron less
      :return: A conceptualDftContainer with local and global cDFT properties
      )delim");

  conceptual_dft_submodule.def("calculate_global", &calculateGlobal, pybind11::arg("energy"),
                               pybind11::arg("energy_plus"), pybind11::arg("energy_minus"), R"delim(
      Calculates a set of global conceptual DFT parameters.

      :param energy: The energy of the structure of interest
      :param energy_plus: The energy for the same geometry with one additional electron
      :param energy_minus: The energy for the same geometry with one electron less
      :return: A globalConceptualDftContainer with global cDFT properties
      )delim");

  conceptual_dft_submodule.def("calculate_local", &calculateLocal, pybind11::arg("atomic_charges"),
                               pybind11::arg("atomic_charges_plus"), pybind11::arg("atomic_charges_minus"), R"delim(
      Calculates the condensed to atom Fukui and dual descriptor indices.

      Note: The quality of the resulting indices heavily depends on the quality of the atomic charges.
      We recommend Hirshfeld charges.

      :param: atomic_charges The atomic charges of the structure of interest
      :param: atomic_charges_plus The atomic charges for the same geometry with one additional electron
      :param: atomic_charges_minus The atomic charges for the same geometry with one electron less
      :return: A localConceptualDftContainer with global cDFT properties
      )delim");

  conceptual_dft_submodule.def("calculate_chemical_potential", &calculateChemicalPotential, pybind11::arg("energy"),
                               pybind11::arg("energy_plus"), pybind11::arg("energy_minus"), R"delim(
      Calculates the conceptual DFT chemical potential.

      Parr, R. G.; Pearson, R. G., J. Am. Chem. Soc. 1983, 105 (26), 7512–7516. https://doi.org/10.1021/ja00364a005.

      :param energy: The energy of the structure of interest
      :param energy_plus: The energy for the same geometry with one additional electron
      :param energy_minus: The energy for the same geometry with one electron less
      :return: The chemical potential
      )delim");

  conceptual_dft_submodule.def("calculate_electronegativity", &calculateElectronegativity, pybind11::arg("energy"),
                               pybind11::arg("energy_plus"), pybind11::arg("energy_minus"), R"delim(
      Calculates the Mulliken electronegativity.

      Mulliken, R. S.; J. Chem. Phys. 1934, 2 (11), 782–793. https://doi.org/10.1063/1.1749394.

      :param energy: The energy of the structure of interest
      :param energy_plus: The energy for the same geometry with one additional electron
      :param energy_minus: The energy for the same geometry with one electron less
      :return: The electronegativity
      )delim");

  conceptual_dft_submodule.def("calculate_hardness", &calculateHardness, pybind11::arg("energy"),
                               pybind11::arg("energy_plus"), pybind11::arg("energy_minus"), R"delim(
      Calculates the chemical hardness.

      Parr, R. G.; Pearson, R. G.; J. Am. Chem. Soc. 1983, 105 (26), 7512–7516. https://doi.org/10.1021/ja00364a005.

      :param energy: The energy of the structure of interest
      :param energy_plus: The energy for the same geometry with one additional electron
      :param energy_minus: The energy for the same geometry with one electron less
      :return: The chemical hardness.
      )delim");

  conceptual_dft_submodule.def("calculate_softness", &calculateSoftness, pybind11::arg("energy"),
                               pybind11::arg("energy_plus"), pybind11::arg("energy_minus"), R"delim(
      Calculates the chemical softness.

      Yang, W.; Parr, R. G.; PNAS 1985, 82 (20), 6723–6726. https://doi.org/10.1073/pnas.82.20.6723.

      :param energy: The energy of the structure of interest
      :param energy_plus: The energy for the same geometry with one additional electron
      :param energy_minus: The energy for the same geometry with one electron less
      :return: The chemical softness.
      )delim");

  conceptual_dft_submodule.def("calculate_electrophilicity", &calculateElectrophilicity, pybind11::arg("energy"),
                               pybind11::arg("energy_plus"), pybind11::arg("energy_minus"), R"delim(
      Calculates the electrophilicity.

      Parr, R. G.; Szentpály, L. v.; Liu, S.; J. Am. Chem. Soc. 1999, 121 (9), 1922–1924.
      https://doi.org/10.1021/ja983494x.

      :param energy: The energy of the structure of interest
      :param energy_plus: The energy for the same geometry with one additional electron
      :param energy_minus: The energy for the same geometry with one electron less
      :return: The electrophilicity.
      )delim");

  conceptual_dft_submodule.def("calculate_fukui_plus", &calculateFukuiPlus, pybind11::arg("atomic_charges"),
                               pybind11::arg("atomic_charges_plus"), pybind11::arg("atomic_charges_minus"), R"delim(
      Calculates the Fukui indices for nucleophilic attacks.
      Note: The quality of the results heavily depends on the quality of the atomic charges.
      We recommend Hirshfeld charges.

      Parr, R. G.; Yang, W.; J. Am. Chem. Soc. 1984, 106 (14), 4049–4050. https://doi.org/10.1021/ja00326a036.

      :param: atomic_charges The atomic charges of the structure of interest
      :param: atomic_charges_plus The atomic charges for the same geometry with one additional electron
      :param: atomic_charges_minus The atomic charges for the same geometry with one electron less
      :return: The Fukui indices for nucleophilic attacks.
      )delim");

  conceptual_dft_submodule.def("calculate_fukui_minus", &calculateFukuiMinus, pybind11::arg("atomic_charges"),
                               pybind11::arg("atomic_charges_plus"), pybind11::arg("atomic_charges_minus"), R"delim(
      Calculates the Fukui indices for nucleophilic attacks.
      Note: The quality of the results heavily depends on the quality of the atomic charges.
      We recommend Hirshfeld charges.

      Parr, R. G.; Yang, W.; J. Am. Chem. Soc. 1984, 106 (14), 4049–4050. https://doi.org/10.1021/ja00326a036.

      :param: atomic_charges The atomic charges of the structure of interest
      :param: atomic_charges_plus The atomic charges for the same geometry with one additional electron
      :param: atomic_charges_minus The atomic charges for the same geometry with one electron less
      :return: The Fukui indices for electrophilic attacks.
      )delim");

  conceptual_dft_submodule.def("calculate_fukui_radical", &calculateFukuiRadical, pybind11::arg("atomic_charges"),
                               pybind11::arg("atomic_charges_plus"), pybind11::arg("atomic_charges_minus"), R"delim(
      Calculates the Fukui indices for radical attacks.
      Note: The quality of the results heavily depends on the quality of the atomic charges.
      We recommend Hirshfeld charges.
      The relevance of the radical Fukui function is known to be limited. It is included here for completeness only.

      Parr, R. G.; Yang, W.; J. Am. Chem. Soc. 1984, 106 (14), 4049–4050. https://doi.org/10.1021/ja00326a036.

      :param: atomic_charges The atomic charges of the structure of interest
      :param: atomic_charges_plus The atomic charges for the same geometry with one additional electron
      :param: atomic_charges_minus The atomic charges for the same geometry with one electron less
      :return: The Fukui indices for radical attacks.
      )delim");

  conceptual_dft_submodule.def("calculate_dual_descriptor", &calculateDualDescriptor, pybind11::arg("atomic_charges"),
                               pybind11::arg("atomic_charges_plus"), pybind11::arg("atomic_charges_minus"), R"delim(
      Calculates the the condensed to atom dual descriptor.
      Note: The quality of the results heavily depends on the quality of the atomic charges.
      We recommend Hirshfeld charges.

      Morell, C.; Grand, A.; Toro-Labbé, A.; J. Phys. Chem. A 2005, 109 (1), 205–212. https://doi.org/10.1021/jp046577a.

      :param: atomic_charges The atomic charges of the structure of interest
      :param: atomic_charges_plus The atomic charges for the same geometry with one additional electron
      :param: atomic_charges_minus The atomic charges for the same geometry with one electron less
      :return: The condensed-to-atom dual descriptor.
      )delim");
}
