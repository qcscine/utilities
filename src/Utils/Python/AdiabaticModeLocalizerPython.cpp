/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/GeometricDerivatives/AdiabaticModeLocalizer.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <boost/optional.hpp>

using namespace Scine::Utils;

void init_adiabatic_mode_localizer(pybind11::module& m) {
  pybind11::class_<AdiabaticModesContainer> adiabaticModesContainer(m, "AdiabaticModesContainer", R"delim(
    A class storing the localized adiabatic properties.
    This class is an extension of the standard normal mode container, with modes being accessed by specifying the
    internal coordinate of interest.
    Note: Currently only interatomic stretching coordinates, e.g., bonds, are supported as internal coordinates.
    )delim");

  adiabaticModesContainer.def("add_mode", &AdiabaticModesContainer::addMode, pybind11::arg("internal_coordinate"),
                              pybind11::arg("mode"), pybind11::arg("force_constant"), R"delim(
    Adds a mode to the container.

    :param internal_coordinate: The internal coordinate to which this mode corresponds.
    :param mode: The mode to be added.
    :param force_constant: The adiabatic force constant of the mode.
    )delim");

  adiabaticModesContainer.def("get_mode", &AdiabaticModesContainer::getMode, pybind11::arg("internal_coordinate"), R"delim(
    Returns the localized vibrational mode corresponding to the given internal coordinate.

    :param: internal_coordinate The internal coordinate of interest.
    :return: The localized mode.
    )delim");

  adiabaticModesContainer.def("get_mode_as_molecular_trajectory", &AdiabaticModesContainer::getModeAsMolecularTrajectory,
                              pybind11::arg("internal_coordinate"), pybind11::arg("atoms"),
                              pybind11::arg("scaling_factor"), R"delim(
    Returns a molecular trajectory corresponding to a vibrational mode.

    :param internal_coordinate: The internal coordinate of interest.
    :param atoms: The atom collection of interest.
    :param scaling_factor: The scaling factor applied to the mode to obtain the maximum displacement.
    :return: The molecular trajectory representing the mode.
    )delim");

  adiabaticModesContainer.def("get_wave_numbers", &AdiabaticModesContainer::getWaveNumbers, R"delim(
    Gets the wave numbers corresponding to the localized modes.
    
    :return: Mapping between internal coordinates and adiabatic wave numbers.
    )delim");

  adiabaticModesContainer.def("get_force_constants", &AdiabaticModesContainer::getForceConstants, R"delim(
    Gets the force constants corresponding to the localized modes.
    
    :return: Mapping between internal coordinates and adiabatic force constants.
    )delim");

  pybind11::class_<AdiabaticModeLocalizer> adiabaticModeLocalizer(m, "AdiabaticModeLocalizer", R"delim(
    This class calculates the localized vibrational modes and adiabatic force constants corresponding to internal
    coordinates based on the local vibrational mode theory by Konkoli and Cremer. Adiabatic force constants are
    equivalent to the relaxed force constants derived from Decius' compliance matrix.
  
    Konkoli, Z.; Cremer, D. Int. J. Quantum Chem. 1998, 67 (1), 1–9.
    https://doi.org/10.1002/(SICI)1097-461X(1998)67:1<1::AID-QUA1>3.0.CO;2-Z.
  
    Implemented after:
    Kraka, E.; Zou, W.; Tao, Y. Wiley Interdiscip. Rev. Comput. Mol. Sci. 2020, e1480. https://doi.org/10.1002/wcms.1480.

    Note: Currently only interatomic stretching coordinates, e.g., bonds, are supported as internal coordinates.
    )delim");

  adiabaticModeLocalizer.def(pybind11::init<Eigen::MatrixXd, AtomCollection, BondOrderCollection, double>(),
                             pybind11::arg("hessian"), pybind11::arg("atoms"), pybind11::arg("bond_orders"),
                             pybind11::arg("bonding_threshold") = 0.5, R"delim(
    Constructor for the AdiabaticModeLocalizer class processing a bond order collection.
    All bonds with a bond order larger than bonding_threshold will be used as internal coordinates.
   
    :param hessian: The cartesian, not mass-weighted Hessian matrix.
    :param atoms: The atom collection of interest.
    :param bond_orders: The bond order collection of the structure of interest.
    :param: bonding_threshold The bond order threshold beyond which bonds are analyzed in terms of localized modes. (default: 0.5)
    )delim");

  adiabaticModeLocalizer.def(pybind11::init<Eigen::MatrixXd, AtomCollection, const std::vector<std::pair<int, int>>>(),
                             pybind11::arg("hessian"), pybind11::arg("atoms"), pybind11::arg("internal_coordinates"), R"delim(
    Constructor for the AdiabaticModeLocalizer class processing indices of atom pairs as internal coordinates.
    
    :param hessian: The cartesian, not mass-weighted Hessian matrix.
    :param atoms: The atom collection of interest.
    :param internal_coordinates:  The internal coordinates of interest.
    )delim");

  adiabaticModeLocalizer.def("localize_modes", &AdiabaticModeLocalizer::localizeModes, R"delim(
    Calculates the adiabatic localized modes corresponding to the given internal coordinates.
   
    References:
    Konkoli, Z.; Cremer, D. Int. J. Quantum Chem. 1998, 67 (1), 1–9.
    https://doi.org/10.1002/(SICI)1097-461X(1998)67:1<1::AID-QUA1>3.0.CO;2-Z.
   
    Kraka, E.; Zou, W.; Tao, Y. Wiley Interdiscip. Rev. Comput. Mol. Sci. 2020, e1480.
    https://doi.org/10.1002/wcms.1480.
   
    :return: A container with the localized adiabatic modes, force constants and wavenumbers.
    )delim");
}
