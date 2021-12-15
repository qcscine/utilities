/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_PERIODICBOUNDARIES_H
#define UTILS_PERIODICBOUNDARIES_H

#include "Utils/Bonds/BondOrderCollection.h"
#include "Utils/Constants.h"
#include "Utils/Typenames.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>

namespace Scine {
namespace Utils {

/*!
 * Class defining Periodic Boundaries
 * The class can be initialized either directly via a matrix defining the imaged cell or via the lengths and angles of
 * the unit cell.
 * The cell matrix is defined with following properties:
 * vector a is defined to be along the x-axis
 * vector b is defined to be in the x-y-plane
 * The transformation matrix is defined as
 * (a0 a1 a2
 *  b0 b1 b2
 *  c0 c1 c2)
 *
 *  alpha = <) b and c
 *  beta = <) a and c
 *  gamma = <) a and b
 *
 * The class can then transform Cartesian coordinates into Relative coordinates and vice versa and positions can be
 * imaged into the cell.
 */
class PeriodicBoundaries {
 public:
  explicit PeriodicBoundaries(double cubeLength = 1.0);
  PeriodicBoundaries(const PeriodicBoundaries& rhs);
  explicit PeriodicBoundaries(Eigen::Matrix3d matrix);

  explicit PeriodicBoundaries(const Eigen::Vector3d& lengths, const Eigen::Vector3d& angles, bool isBohr = true,
                              bool isDegrees = true);

  explicit PeriodicBoundaries(std::string periodicBoundariesString, const std::string& delimiter = ",",
                              bool isBohr = true, bool isDegrees = true);

  /* Transformation Block */
  /**
   * @brief Transform between Cartesian and Relative Coordinates and returns new PositionCollection as copy
   * @param positions A PositionCollection
   * @param relativeToCartesian Whether Relative Coordinates are transformed into Cartesian or vice versa
   * @return PositionCollection Transformed PositionCollection
   */
  PositionCollection transform(const PositionCollection& positions, bool relativeToCartesian = true) const;

  /**
   * @brief Transform between Cartesian and Relative Coordinates and returns new Position as copy
   * @param position A single Position
   * @param relativeToCartesian Whether Relative Coordinates are transformed into Cartesian or vice versa
   * @return Position Transformed Position
   */
  Position transform(const Position& position, bool relativeToCartesian = true) const;

  /**
   * @brief Transform between Cartesian and Relative Coordinates in place
   * @param positions A PositionCollection
   * @param relativeToCartesian Whether Relative Coordinates are transformed into Cartesian or vice versa
   * @return void
   */
  inline void transformInPlace(PositionCollection& positions, bool relativeToCartesian = true) const {
    if (relativeToCartesian) {
      positions *= _matrix;
    }
    else {
      positions *= _invMatrix;
    }
  };

  /**
   * @brief Transform between Cartesian and Relative Coordinates in place
   * @param position A single Position
   * @param relativeToCartesian Whether Relative Coordinates are transformed into Cartesian or vice versa
   * @return void
   */
  inline void transformInPlace(Eigen::Ref<Position> position, bool relativeToCartesian = true) const {
    if (position.cols() != 3) {
      throw std::runtime_error("The given position(s) are not of dimension 3. This is not possible.");
    }
    if (relativeToCartesian) {
      position *= _matrix;
    }
    else {
      position *= _invMatrix;
    }
  }

  /* Translation Block */
  /**
   * @brief Translate single Position in Cartesian Coordinates into unit cell and returns new Position as copy.
   * @param position Position to be translated
   * @param relShift An optional shift vector in Relative coordinates to apply after the translation into the cell
   * @return Position Translated Position
   */
  Position translatePositionsIntoCell(const Position& position,
                                      const Eigen::RowVector3d& relShift = Eigen::RowVector3d::Zero()) const;

  /**
   * @brief Translate PositionCollection in Cartesian Coordinates into unit cell and returns new Collection as copy.
   * @param positions PositionCollection to be translated
   * @param relShift An optional shift vector in Relative coordinates to apply after the translation into the cell
   * @return PositionCollection Translated PositionCollection
   */
  PositionCollection translatePositionsIntoCell(const PositionCollection& positions,
                                                const Eigen::RowVector3d& relShift = Eigen::RowVector3d::Zero()) const;

  /**
   * @brief Translate single Position in Cartesian Coordinates into unit cell in place
   * @param position Position to be translated
   * @param relShift An optional shift vector in Relative coordinates to apply after the translation into the cell
   * @return void
   */
  void translatePositionsIntoCellInPlace(Eigen::Ref<Position> position,
                                         const Eigen::RowVector3d& relShift = Eigen::RowVector3d::Zero()) const;

  /**
   * @brief Translate PositionCollection in Cartesian Coordinates into unit cell in place
   * @param positions PositionCollection to be translated
   * @param relShift An optional shift vector in Relative coordinates to apply after the translation into the cell
   * @return void
   */
  void translatePositionsIntoCellInPlace(PositionCollection& positions,
                                         const Eigen::RowVector3d& relShift = Eigen::RowVector3d::Zero()) const;

  /**
   * @brief Get the squared distance between two Positions using the minimum image convention with a faster algorithm,
   * which is only valid if the distance within the cell is less than half of the shortest distance between two opposing
   * faces of the cell
   *
   * @note use getDistanceSquared function in GeometryUtilities to automatically decide whether this is reasonable.
   *
   * @param p1 First Position
   * @param p2 Second Position
   * @return double squared Distance
   */
  double fastMinimumImageDistanceSquared(const Position& p1, const Position& p2) const;

  /**
   * @brief Get the squared distance between two Positions using the minimum image convention calculating all image
   * distances
   *
   * @note use getDistanceSquared function in GeometryUtilities to automatically decide whether this is necessary.
   *
   * @param p1 First Position
   * @param p2 Second Position
   * @return double squared Distance
   */
  double bruteForceMinimumImageDistanceSquared(const Position& p1, const Position& p2) const;

  /**
   * @brief Get the vector of minimum length between two Positions based on a brute-force image distances
   *
   * @param p1 First Position
   * @param p2 Second Position
   * @return The displacement vector
   *
   * @note The vector is only guaranteed to be unique if its length is smaller than half of the shortest distance
   * between two opposing faces of the cell.
   */
  Displacement bruteForceMinimumImageDisplacementVector(const Position& p1, const Position& p2) const;

  /**
   * @brief Get the squared distances between a Position and all images of another Position.
   * @param p1 First Position
   * @param p2 Second Position
   * @return std::vector<double> squared Distances
   */
  std::vector<double> getAllImageDistancesSquared(const Position& p1, Eigen::Ref<Position> p2) const;

  /**
   * @brief Gets the displacement vectors pointing from a Position to all images of another Position.
   * @param p1 First Position
   * @param p2 Second Position
   * @return std::vector<double> displacement vectors
   */
  std::vector<Displacement> getAllImageDisplacementVectors(const Position& p1, Eigen::Ref<Position> p2) const;

  /**
   * @brief Determine whether the smallest distance between two positions is between their real space positions or
   * via an image.
   *
   * @return true Closest distance via image.
   * @return false Closest distance between position within the periodic boundaries.
   */
  bool minimumDistanceViaImage(Position p1, Position p2) const;

  /**
   * @brief Whether a given position is within the periodic boundaries.
   *
   * @param pos The position.
   * @return bool
   */
  bool isWithinCell(const Position& position) const;
  /**
   * @brief Whether all positions in a PositionCollection are within the periodic boundaries.
   *
   * @param pos The positions.
   * @return bool
   */
  bool isWithinCell(const PositionCollection& positions) const;

  inline std::string getPeriodicBoundariesString() const {
    return std::to_string(_aNorm) + "," + std::to_string(_bNorm) + "," + std::to_string(_cNorm) + "," +
           std::to_string(_alpha) + "," + std::to_string(_beta) + "," + std::to_string(_gamma);
  };

  inline Eigen::Matrix3d getCellMatrix() const {
    return _matrix;
  };

  inline void setCellMatrix(Eigen::Matrix3d matrix) {
    _matrix = std::move(matrix);
    setMembers();
  };

  inline Eigen::Vector3d getA() const {
    return _matrix.row(0);
  };

  inline Eigen::Vector3d getB() const {
    return _matrix.row(1);
  };

  inline Eigen::Vector3d getC() const {
    return _matrix.row(2);
  };

  inline double getAlpha() const {
    return _alpha;
  };

  inline double getBeta() const {
    return _beta;
  };

  inline double getGamma() const {
    return _gamma;
  };

  inline double getSmallestPerpendicularSquared() const {
    return _minPerpendicularSq;
  };

  inline double getBiggestPossibleDistanceSquared() const {
    return _biggestDistanceSquared;
  };

  PeriodicBoundaries& operator=(const PeriodicBoundaries& rhs);
  PeriodicBoundaries& operator=(const Eigen::Matrix3d& rhs);
  /**
   * @brief Operator overload, increase periodic boundaries
   * @param scalingFactors The scaling factors for a, b, and c vectors
   * @return The bigger periodic boundaries
   */
  PeriodicBoundaries operator*(const Eigen::Vector3d& scalingFactors) const;
  /**
   * @brief Operator overload, increase periodic boundaries in place
   * @param scalingFactors The scaling factors for a, b, and c vectors
   * @return The bigger periodic boundaries
   */
  PeriodicBoundaries& operator*=(const Eigen::Vector3d& scalingFactors);
  /**
   * @brief Operator overload, increase periodic boundaries
   * @param scalingFactor The scaling factor for all 3 dimensions
   * @return The bigger periodic boundaries
   */
  PeriodicBoundaries operator*(double scalingFactor) const;
  /**
   * @brief Operator overload, increase periodic boundaries in place
   * @param scalingFactor The scaling factors for all 3 dimensions
   * @return The bigger periodic boundaries
   */
  PeriodicBoundaries& operator*=(double scalingFactor);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW // See http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html

      private : Eigen::Matrix3d _matrix;
  Eigen::Matrix3d _invMatrix;
  const double _eps = 1e-6;
  /* lengths */
  double _aNorm;
  double _bNorm;
  double _cNorm;
  /* angles */
  double _alpha;
  double _beta;
  double _gamma;
  /* Cell distance properties */
  double _biggestDistanceSquared;
  double _minPerpendicularSq;

  // @brief Convenience functions to set close-to-zero-values to zero
  void reduceNoise(Eigen::Vector3d& vec) const;

  void reduceNoise(Eigen::Matrix3d& mat) const;

  void setMembers();

  // @brief Convenience functions for angle conversions
  static inline double getRadians(const double angle) {
    return angle * Constants::pi / 180.0;
  };

  static inline double getDegrees(const double angle) {
    return angle * 180.0 / Constants::pi;
  };

  void constructMembersFromLengthsAndAngles(const Eigen::Vector3d& lengths, const Eigen::Vector3d& angles, bool isBohr,
                                            bool isDegrees);
};

bool operator==(const PeriodicBoundaries& lhs, const PeriodicBoundaries& rhs);
bool operator!=(const PeriodicBoundaries& lhs, const PeriodicBoundaries& rhs);

} // namespace Utils
} // namespace Scine
#endif // UTILS_PERIODICBOUNDARIES_H
