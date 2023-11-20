/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_SPGINTERFACE_H_
#define UTILS_SPGINTERFACE_H_

#include "Utils/Geometry/PeriodicSystem.h"
#include "Utils/Typenames.h"
#include <Eigen/Core>
#include <functional>
#include <memory>
#include <tuple>
#include <unordered_set>
#include <utility>

namespace Scine {
namespace Utils {

/**
 * @class SpgInterface SpgInterface.h
 * @brief A class to ease calls to the SPGLib C-API
 */
class SpgInterface {
 public:
  /// @brief Static functions only.
  SpgInterface() = delete;
  /**
   * @brief Get a PeriodicSystem that represents the primitive cell of the given system.
   * @param system The system we want the primitive cell of.
   * @param epsilon The symmetry and equal position precision.
   * @param solidStateOnly Whether you want only a primitive cell of all the solid state atoms, by default false.
   * @return PeriodicSystem The system containing only the primitive cell.
   */
  static PeriodicSystem findPrimitiveCellSystem(const PeriodicSystem& system, double epsilon, bool solidStateOnly = false);
  /**
   * @brief Get a list of matrices that are symmetry equivalent to the given matrix.
   * @param cellMatrix The cell matrix we want equivalents of.
   * @param epsilon The symmetry and equal position precision.
   * @return std::vector<Eigen::Matrix3d> The list of matrices.
   */
  static std::vector<Eigen::Matrix3d> findAlternativePbcs(const Eigen::Matrix3d& cellMatrix, double epsilon);

  /**
   * @struct SpgInterface::SymmetryOperation SpgInterface.h
   * @brief A struct representing a symmetry operations (rotation + translation in fractional coordinates)
   */
  struct SymmetryOperation {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // See http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
        Eigen::Matrix3d rotation = Eigen::Matrix3d::Identity(); // transform to double for easier multiplication
    Eigen::RowVector3d translation = Eigen::RowVector3d ::Zero();
    Eigen::Matrix3d pbc;
    inline SymmetryOperation(int cRotation[3][3], double cTranslation[3], const Eigen::Matrix3d& pbc) {
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          rotation(j, i) = static_cast<double>(cRotation[i][j]); // shifted indices are on purpose
        }
        translation[i] = cTranslation[i];
      }
      this->pbc = pbc;
    };
  };

  /**
   * @struct SpgInterface::CppCell SpgInterface.h
   * @brief A class holding information required by SpgLib in more convenient C++ data structures + basic utilities.
   */
  class CppCell {
   public:
    PeriodicBoundaries lattice;
    PositionCollection positions;
    std::vector<int> types;
    inline CppCell(const PeriodicBoundaries& a, PositionCollection b, std::vector<int> c)
      : lattice(a), positions(std::move(b)), types(std::move(c)){};
    /**
     * @brief Whether the two cells are approximately equal including various symmetry operations.
     * @param rhs A different cell.
     * @param epsilon The symmetry and equal position precision.
     * @return bool If equal or not
     */
    bool isApprox(CppCell rhs, double epsilon) const;
    /**
     * @brief Whether the two cells are approximately equal with predefined symmetry operations.
     * @param rhs A different cell.
     * @param epsilon The symmetry and equal position precision.
     * @param lhsOperations The List of operations on the executing class.
     * @param rhsOperations The List of operations on the compared class.
     * @return bool If equal or not
     */
    bool isApprox(CppCell rhs, double epsilon, const std::vector<SymmetryOperation>& lhsOperations,
                  const std::vector<SymmetryOperation>& rhsOperations) const;

   private:
    bool isApproxImpl(const Scine::Utils::SpgInterface::CppCell& rhs, double epsilon,
                      const std::vector<PositionCollection>& lhsPositions,
                      const std::vector<PositionCollection>& rhsPositions) const;
  };
  /**
   * @brief Get a PeriodicSystem that represents the primitive cell of the given system.
   * @param system The system we want the primitive cell of.
   * @param epsilon The symmetry and equal position precision.
   * @param solidStateOnly Whether you want only a primitive cell of all the solid state atoms, by default false.
   * @return PeriodicSystem The system containing only the primitive cell.
   */

  /**
   * @brief Get a map of indices of positions in lhs to all positions in rhs
   * @param lhs A cell we want the assignment from.
   * @param rhs A cell we want the assignments to.
   * @param epsilon The symmetry and equal position precision.
   * @return std::map<int, std::vector<int>> The map from left to right
   */
  static std::map<int, std::vector<int>> mapSymmetryEquivalentPositions(const CppCell& lhs, const CppCell& rhs, double epsilon);
  /**
   * @brief Get a Cell that represents the primitive cell of the given system.
   * @param system The system we want the primitive cell of.
   * @param epsilon The symmetry and equal position precision.
   * @param solidStateOnly Whether you want only a primitive cell of all the solid state atoms, by default false.
   * @return CppCell The cell containing only the primitive cell.
   */
  static CppCell findPrimitiveCell(const PeriodicSystem& system, double epsilon, bool solidStateOnly = false);
  /**
   * @brief Get a list of symmetry operations for the given cell.
   * @param cell The cell we want the symmetry operations from.
   * @param epsilon The symmetry and equal position precision.
   * @return std::vector<SymmetryOperation> The list of operations.
   */
  static std::vector<SymmetryOperation> findSymmetryOperations(const CppCell& cell, double epsilon);
  /**
   * @brief Get a list of positions that are symmetrically equivalent to the given positions.
   * @param positions The positions.
   * @param operations A list of operations that can be applied to the positions.
   * @return std::vector<PositionCollection> The list of equivalent positions, should contain positions.
   */
  static std::vector<PositionCollection> getSymmetryEquivalents(const PositionCollection& positions,
                                                                const std::vector<SpgInterface::SymmetryOperation>& operations);
  /**
   * @brief Get a list of positions that are symmetrically equivalent to the positions in the given cell.
   * @param cell The Cell specifying positions and periodic boundaries.
   * @param epsilon The symmetry and equal position precision.
   * @return std::vector<PositionCollection> The list of equivalent positions, should contain positions.
   */
  static std::vector<PositionCollection> getSymmetryEquivalents(const CppCell& cell, double epsilon);

 private:
  /**
   * @struct SpgInterface::CCell SpgInterface.h
   * @brief A class holding information required by SpgLib in C data structures (partially wrapped in smart pointers)
   */
  struct CCell {
    double lattice[3][3];
    std::shared_ptr<double (*)[3]> positions;
    std::shared_ptr<int*> types;
    int nAtoms;
    inline CCell(double a[3][3], std::shared_ptr<double (*)[3]> b, std::shared_ptr<int*> c, int d)
      : positions(std::move(b)), types(std::move(c)), nAtoms(d) {
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          lattice[i][j] = a[i][j];
        }
      }
    };
  };

  /// changes positions according to given symmetry operation
  static void applySymmetryOperation(PositionCollection& positions, const SymmetryOperation& operation);
  /// find symmetry, but with CCell input, constructing cell matrix information
  static std::vector<SymmetryOperation> findSymmetryOperations(const CCell& cell, double epsilon);
  /// find symmetry, but with CCell input + cell matrix
  static std::vector<SymmetryOperation> findSymmetryOperations(const CCell& cell, double epsilon, const Eigen::Matrix3d& pbc);
  /// convert a PeriodicSystem to a CCell
  static CCell systemToCell(const PeriodicSystem& system, bool solidStateOnly = false);
  /// convert a CCell to a CppCell
  static CppCell cellToCppCell(const CCell& cell);
  /// convert a CppCell to a CCell
  static CCell cppCellToCell(const CppCell& cell);
  /// convert a PeriodicSytem to a CppCell
  static CppCell systemToCppCell(const PeriodicSystem& system, bool solidStateOnly = false);
  /**
   * Find minimum distance between given position lhsPos within rhsPosition of the same type,
   * taking pbc into consideration, but no position manipulations are carried out.
   */
  static std::pair<int, int> minDistanceAndIndex(int lhsType, const Position& lhsPos, std::vector<int> rhsTypes,
                                                 const PositionCollection& rhsPositions, const PeriodicBoundaries& pbc);
};

} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_SPGINTERFACE_H_
