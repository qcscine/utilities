/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_BSPLINES_REACTIONPROFILEINTERPOLATION_H
#define UTILS_BSPLINES_REACTIONPROFILEINTERPOLATION_H

#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Typenames.h"
#include <Eigen/Dense>
#include <memory>

namespace Scine {
namespace Utils {
namespace BSplines {

/**
 * @brief Data class for one molecular spline, with evaluation functionality.
 */
class TrajectorySpline {
 public:
  /**
   * @brief Construct a new TrajectorySpline object.
   * @param elements The elements of the atoms in the trajectory.
   * @param knots The list of all knot positions in the interval [0.0, 1.0].
   * @param data The data at the knot position dimensions (nKnots, nAtoms*3 + 1)
   *             The energy is store in the 0th index before all structure data.
   * @param tsPosition The position of a transition state in the spline. A valid
   *                   transition state has to be in the intervall [0.0, 1.0].
   *                   Defaults to no transition state (value -1).
   */
  TrajectorySpline(const ElementTypeCollection& elements, const Eigen::VectorXd& knots, const Eigen::MatrixXd& data,
                   double tsPosition = -1.0);
  ~TrajectorySpline() = default;
  /**
   * @brief Evaluate the spline at the given position in the interval [0.0, 1.0]
   * @param position The position to interpolate at.
   * @param degree The degree of the polynomial fit to use in the interpolation.
   * @return std::tuple<double, AtomCollection> The energy and structure at the
   *         requested position.
   */
  std::tuple<double, AtomCollection> evaluate(const double& position, const unsigned int& degree) const;

  const ElementTypeCollection elements;
  const Eigen::VectorXd knots;
  const Eigen::MatrixXd data;
  const double tsPosition;
};

/**
 * @brief Factory for TrajectorySplines of reaction paths.
 *
 * Including interpolation of an energy associated with the molecular structures.
 */
class ReactionProfileInterpolation {
 public:
  ReactionProfileInterpolation() = default;
  ~ReactionProfileInterpolation() = default;
  /**
   * @brief Adds a single structure to the back of the stored trajectory.
   * @param atoms The collection of atoms of the structure to add.
   * @param energy The energy of the structure to add.
   * @param isTS Set to true if the given structure is the transition state
   *             of the trajectory to fit.
   */
  void appendStructure(const AtomCollection& atoms, const double& energy, bool isTS = false);
  /**
   * @brief Deletes all stored data.
   */
  void clear();
  /**
   * @brief Generates a spline using all currently stored data.
   *
   * The first structure in the stored data will correspond to the position 0.0
   * of the spline, the last one to position 1.0.
  _*
   * @param nInterpolationPoints The number of interpolation points used in the
   *                             final spline, values smaller than the number of
   *                             points in the trajectory will result in
   *                             compression and possibly smoothing of the
   *                             curve. The number of points needs to be bigger
   *                             than 3 if a transition state has been given
   *                             or greater than 2 if not.
   * @param degree The degree of fit to use in the compression of the given
   *               trajectory down to the requested number of interpolation
   *               points
   * @return TrajectorySpline A spline object fitting the given trajectory
   */
  TrajectorySpline spline(unsigned int nInterpolationPoints = 11, unsigned int degree = 3) const;
  /**
   * @brief Get the current TS position.
   * @return double The position of the TS in the interpolation interval of
   *                [0.0, 1.0], based on the current data.
   *                Will throw an exception if no TS has been given
   */
  double getCurrentTSPosition() const;
  /// @brief If true will (massweighted) fit each structure to the previous one.
  bool useQuaternionFit = true;

 private:
  std::unique_ptr<AtomCollection> _start = nullptr;
  std::unique_ptr<unsigned int> _tsIdx = nullptr;
  Eigen::MatrixXd _data = Eigen::MatrixXd::Zero(0, 0);
  unsigned int _nPoints = 0;
};

} // namespace BSplines
} // namespace Utils
} // namespace Scine

#endif // UTILS_BSPLINES_REACTIONPROFILEINTERPOLATION_H
