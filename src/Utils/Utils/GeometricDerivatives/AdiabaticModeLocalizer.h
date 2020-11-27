/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_ADIABATICMODELOCALIZER_H
#define UTILS_ADIABATICMODELOCALIZER_H

#include "Utils/Bonds/BondOrderCollection.h"
#include "Utils/Constants.h"
#include "Utils/GeometricDerivatives/NormalModesContainer.h"
#include "Utils/Geometry/AtomCollection.h"
#include <Eigen/Core>
#include <exception>

namespace Scine {
namespace Utils {

namespace AdiabaticLocalizationExceptions {
class modeReadded : public std::exception {
  const char* what() const noexcept final {
    return "There already is a mode corresponding to the given internal coordinate.";
  }
};

class InvalidModeRequested : public std::exception {
  const char* what() const noexcept final {
    return "No mode corresponding to given internal coordinate";
  }
};

class InvalidHessianDimensions : public std::exception {
  const char* what() const noexcept final {
    return "Hessian matrix with incorrect dimensions was passed to the adiabatic mode localizer";
  }
};

} // namespace AdiabaticLocalizationExceptions

/**
 * @class AdiabaticModeContainer @file AdiabaticModeLocalizer.h
 * @brief A class storing the localized adiabatic properties
 *
 * This class is an extension of the standard normal mode container, with modes being accessed by specifying the
 * internal coordinate of interest.
 *
 * @note Currently only interatomic stretching coordinates, e.g., bonds, are supported as internal coordinates.
 */
class AdiabaticModesContainer {
 public:
  /**
   * @brief Adds a mode to the container.
   * @param internalCoordinate The internal coordinate to which this mode corresponds.
   * @param mode The mode to be added.
   * @param forceConstant The adiabatic force constant of the mode.
   */
  void addMode(std::pair<int, int> internalCoordinate, NormalMode mode, double forceConstant);

  /**
   * @brief Returns a const reference of the localized vibrational mode corresponding to the given internalCoordinate.
   * @param internalCoordinate The internal coordinate of interest.
   * @return The localized mode.
   */
  const DisplacementCollection& getMode(std::pair<int, int> internalCoordinate) const;

  /**
   * @brief Returns a molecular trajectory corresponding to a vibrational mode.
   *
   * @param internalCoordinate The internal coordinate of interest.
   * @param atoms The AtomCollection of interest.
   * @param scalingFactor The scaling factor applied to the mode to obtain the maximum displacement.
   * @return MolecularTrajectory The molecular trajectory representing the mode.
   */
  MolecularTrajectory getModeAsMolecularTrajectory(std::pair<int, int> internalCoordinate,
                                                   const Utils::AtomCollection& atoms, double scalingFactor) const;

  /**
   * @brief Gets the wave numbers corresponding to the localized modes.
   * @return std::map<std::pair<int,int>,double> Mapping between internal coordinates and adiabatic wave numbers.
   */
  std::map<std::pair<int, int>, double> getWaveNumbers();

  /**
   * @brief Gets the force constants corresponding to the localized modes.
   * @return std::map<std::pair<int,int>,double> Mapping between internalCoordinates and adiabatic force constants.
   */
  std::map<std::pair<int, int>, double> getForceConstants() const;

 private:
  /// @brief The underlying normal modes container
  NormalModesContainer modesContainer_;
  /// @brief Mapping between atom pair and storage index
  std::map<std::pair<int, int>, int> storageOrder_;
  /// @brief Map between coordinates and force constants;
  std::map<std::pair<int, int>, double> forceConstants_;
  /// @brief Map between coordinates and wavenumbers;
  std::map<std::pair<int, int>, double> wavenumbers_;

  /**
   * @brief Gets the storage index corresponding to a internalCoordinate.
   *
   * Takes care about the index ordering in the internalCoordinate pair and checks whether the information
   * corresponding to the given indices exists.
   *
   * @param internalCoordinate The internal coordinate of interest.
   * @return The storage index corresponding to the stretching movement between the given atoms.
   */
  int getStorageIndex(std::pair<int, int> internalCoordinate) const;

}; // AdiabaticModeContainer

/**
 * @class AdiabaticModeLocalizer @file AdiabaticModeLocalizer.h
 * @brief Calculates the localized vibrational modes and adiabatic force constants corresponding to internal coordinates
 *        of a molecule
 *
 * This class calculates the localized vibrational modes and adiabatic force constants corresponding to internal
 * coordinates based on the local vibrational mode theory by Konkoli and Cremer. Adiabatic force constants are
 * equivalent to the relaxed force constants derived from Decius' compliance matrix.
 *
 * Konkoli, Z.; Cremer, D. Int. J. Quantum Chem. 1998, 67 (1), 1–9.
 * https://doi.org/10.1002/(SICI)1097-461X(1998)67:1<1::AID-QUA1>3.0.CO;2-Z.
 *
 * Implemented after:
 * Kraka, E.; Zou, W.; Tao, Y. Wiley Interdiscip. Rev. Comput. Mol. Sci. 2020, e1480. https://doi.org/10.1002/wcms.1480.
 *
 * @note Currently only interatomic stretching coordinates, e.g., bonds, are supported as internal coordinates.
 */
class AdiabaticModeLocalizer {
 public:
  /**
   * @brief Constructor for the AdiabaticModeLocalizer class processing a BondOrderCollection.
   *
   * All bonds with a bond order larger than bondingThreshold will be used as internal coordinates.
   *
   * @param hessian The cartesian, not mass-weighted Hessian matrix.
   * @param atoms The atom collection of interest.
   * @param bondOrders The bond order collection of the structure of interest.
   * @param bondingThreshold The bond order threshold beyond which bonds are analyzed in terms of localized modes.
   * (default: 0.5)
   */
  AdiabaticModeLocalizer(const Eigen::MatrixXd& hessian, const AtomCollection& atoms,
                         const BondOrderCollection& bondOrders, const double bondingThreshold = 0.5);

  /**
   * @brief Constructor for the AdiabaticModeLocalizer class processing indices of atom pairs as internal coordinates.
   *
   * @param hessian The cartesian, not mass-weighted Hessian matrix.
   * @param atoms The atom collection of interest.
   * @param internalCoordinates  The internal coordinates of interest.
   */
  AdiabaticModeLocalizer(const Eigen::MatrixXd& hessian, const AtomCollection& atoms,
                         const std::vector<std::pair<int, int>>& internalCoordinates);

  /**
   * @brief Calculates the adiabatic localized modes corresponding to the given internal coordinates.
   *
   * References:
   * Konkoli, Z.; Cremer, D. Int. J. Quantum Chem. 1998, 67 (1), 1–9.
   * https://doi.org/10.1002/(SICI)1097-461X(1998)67:1<1::AID-QUA1>3.0.CO;2-Z.
   *
   * Kraka, E.; Zou, W.; Tao, Y. Wiley Interdiscip. Rev. Comput. Mol. Sci. 2020, e1480.
   * https://doi.org/10.1002/wcms.1480.
   *
   * @return A container with the localized adiabatic modes, force constants and wavenumbers.
   */
  AdiabaticModesContainer localizeModes();

 private:
  /**
   * @brief Calculates the Wilson B matrix for the given internalCoordinates
   *
   * @note The B matrix is also calculated in the irc library.
   */
  void calculateStretchingB();

  Utils::AtomCollection atoms_;
  Eigen::MatrixXd hessian_;                         // Cartesian, not mass-weighted hessian
  std::vector<std::pair<int, int>> internalCoords_; // The indices of the bonds of interest
  Eigen::MatrixXd B_;                               // Wilson B matrix for bonds
  Eigen::MatrixXd M_;                               // Diagonal mass matrix
};                                                  // AdiabaticModeLocalizer

} // namespace Utils
} // namespace Scine

#endif // UTILS_ADIABATICMODELOCALIZER_H
