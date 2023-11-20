/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILSOS_BASISSET_H
#define UTILSOS_BASISSET_H

#include <Utils/DataStructures/Shell.h>
#include <Utils/DataStructures/ShellPairs.h>
#include <Utils/Geometry/AtomCollection.h>
#include <memory>
#include <utility>

namespace Scine {
namespace Utils {
namespace Integrals {

/**
 * @class BasisSet
 * @file BasisSet.h
 * @brief Slightly decorated std::vector of Shell objects.
 * @note Should only be used together with the IntegralCalculator.
 * The main functionality is the automatic construction by a basis set name via the Libint interface.
 */
class BasisSet : public std::vector<Shell> {
 private:
  std::string _name;
  AtomCollection _atoms;

 public:
  const AtomCollection& getAtoms() const;

 private:
  std::size_t _id; // id is used for efficient comparison of two BasisSet objects.
  static int _idProvider;
  std::shared_ptr<ShellPairs> _shellPairs;
  bool _shellPairsEvaluated = false;
  bool _isPureSpherical;

 public:
  /**
   * Returns true of shell pairs are evaluated.
   */
  bool areShellPairsEvaluated() const;

  /**
   * @brief returns the ShellPairs object. Must be evaluated via the IntegralCalculator first!!
   */
  const std::shared_ptr<ShellPairs>& getShellPairs() const;
  /**
   * @brief Setter for ShellPairs object.
   */
  void setShellPairs(std::shared_ptr<ShellPairs>& shellPairs);

  /**
   * @brief Setter for ShellPairs object.
   */
  void deleteShellPairs();

  BasisSet() : _id(++_idProvider){};

  /**
   * @brief Constructor only populated name and atoms. Does not construct the basis set with the name!
   *        This must be done via the interface.
   */
  BasisSet(std::string name, AtomCollection atoms)
    : _name(std::move(name)), _atoms(std::move(atoms)), _id(++_idProvider){};

  BasisSet(AtomCollection atoms) : _atoms(std::move(atoms)), _id(++_idProvider){};

  BasisSet& operator=(const BasisSet&) = default;

  BasisSet(const BasisSet& basis) = default;

  ~BasisSet() = default;

  inline bool operator==(const BasisSet& rhs) const {
    return (this->_id == rhs._id);
  }
  inline bool operator!=(const BasisSet& rhs) const {
    return (this->_id != rhs._id);
  }

  /**
   * @brief Returns true if only solid harmonics are used. Else only cartesian Gaussians are used.
   */
  inline auto isPureSpherical() const -> bool {
    return _isPureSpherical;
  }

  /**
   * @brief Returns true if only solid harmonics are used. Else only cartesian Gaussians are used.
   */
  inline auto setPureSpherical(bool isPureSpherical) -> void {
    _isPureSpherical = isPureSpherical;
  }

  /**
   * @brief Returns the number of basis functions.
   */
  inline auto nbf() const -> std::size_t {
    std::size_t nbf = 0;
    for (auto const& shell : *this) {
      nbf += shell.size();
    }
    return nbf;
  }
  /**
   * @brief Returns maximum angular momentum.
   */
  inline auto max_l() const -> std::size_t {
    std::size_t max_l = 0;
    for (auto const& shell : *this) {
      max_l = std::max(shell.l(), max_l);
    }
    return max_l;
  }
  /**
   * @brief Returns maximum number of primitives, i.e., maximum contraction length.
   */
  inline auto max_nprim() const -> std::size_t {
    std::size_t max_nprim = 0;
    for (auto const& shell : *this) {
      max_nprim = std::max(shell.nprim(), max_nprim);
    }
    return max_nprim;
  }

  /**
   * @brief Returns the map from shell index to index of the first basis function from this shell
   */
  inline auto shell2bf() const -> std::vector<std::size_t> {
    std::vector<std::size_t> result;
    result.reserve(this->size());

    std::size_t n = 0;
    for (auto const& shell : *this) {
      result.push_back(n);
      n += shell.size();
    }

    return result;
  }

  /**
   * @brief returns the ID for comparison of two basis sets.
   */
  inline auto getID() const -> std::size_t {
    return _id;
  }

  /**
   * @brief Moves the whole basis set to new position
   */
  inline auto move(const Displacement& newShift) -> void {
    _id = (++_idProvider);
    for (auto& shell : *this) {
      shell.move(shell.getShift() + newShift);
    }
  }

  /**
   * @brief Moves the whole basis set to new position
   */
  inline auto rotate(const Eigen::Matrix3d& R) -> void {
    _id = (++_idProvider);
    for (auto& shell : *this) {
      shell.move((R * shell.getShift().transpose()).transpose());
    }
  }

  /**
   * @brief append a basis set, and update the atoms member.
   * @param rhs
   * @param rhs_atoms
   */
  auto append(const BasisSet& rhs) -> void;

  /**
   * @brief un-contracts a basis set.
   */
  auto uncontract() -> void;

  /**
   * @brief Computes the map from this object's shells to the corresponding atoms in \c atoms. If no atom matches the
   * origin of a shell, it is mapped to -1.
   * @note coordinates must match \em exactly , i.e. shell2atom[k] == l iff atoms[l].x == *this[k].O[0] && atoms[l].y ==
   * *this[k].O[1] &&  atoms[l].z == *this[k].O[2]
   * @return the map from shell index to the atom in the list \c atoms that coincides with its origin;
   */
  vector<long> shellToAtom(const AtomCollection& atoms) const;

  /**
   * @brief Computes the map from \c shells to the corresponding atoms in \c atoms. Coordinates are compared bit-wise,
   * i.e. shell2atom[k] == l iff atoms[l].x == *this[k].O[0] && atoms[l].y == *this[k].O[1] &&  atoms[l].z ==
   * *this[k].O[2]
   * @param throw_if_no_match If true, and no atom matches the origin of a shell, throw a std::logic_error.
   *        Otherwise such shells will be mapped to -1.
   * @return the map from shell index to the atom in the list \c atoms that coincides with its origin;
   * @throw std::logic_error if throw_if_no_match is true and for at least one shell no matching atom is found.
   */
  static vector<long> shellToAtom(const BasisSet& basis, const AtomCollection& atoms, bool throw_if_no_match);

  /**
   * @brief Computes the map from \c atoms to the corresponding shells in this object. Coordinates are compared bit-wise
   * (@sa BasisSet::shell2atom() )
   * @return the map from atom index to the vector of shell indices whose origins conincide with the atom;
   * @note this does not assume that \c shells are ordered in the order of atoms, as does BasisSet
   */
  std::vector<std::vector<long>> atomToShell(const AtomCollection& atoms) const;

  /**
   * @brief Computes the map from \c atoms to the corresponding shells in \c shells. Coordinates are compared bit-wise
   * (@sa BasisSet::shell2atom() )
   * @return the map from atom index to the vector of shell indices whose origins conincide with the atom;
   * @note this does not assume that \c shells are ordered in the order of atoms, as does BasisSet
   */
  static vector<vector<long>> atomToShell(const AtomCollection& atoms, const BasisSet& basis);
};

} // namespace Integrals
} // namespace Utils
} // namespace Scine

#endif // UTILSOS_BASISSET_H
