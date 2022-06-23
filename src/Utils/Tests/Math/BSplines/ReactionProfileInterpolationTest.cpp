/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/Math/QuaternionFit.h"
#include <Utils/IO/MolecularTrajectoryIO.h>
#include <Utils/Math/BSplines/ReactionProfileInterpolation.h>
#include <Utils/MolecularTrajectory.h>
#include <gmock/gmock.h>
#include <Eigen/Core>

using namespace testing;
namespace Scine {
namespace Utils {
using namespace BSplines;
namespace Tests {

class ReactionProfileInterpolationTest : public Test {
 public:
  void SetUp() override {
  }
};

TEST_F(ReactionProfileInterpolationTest, RPIFailure_Appending_ElementMissMatch) {
  AtomCollection first;
  first.push_back(Atom(ElementType::He, Position(0, 0, 0)));
  first.push_back(Atom(ElementType::Xe, Position(1, 1, 1)));

  AtomCollection second;
  second.push_back(Atom(ElementType::He, Position(0, 0, 0)));
  second.push_back(Atom(ElementType::H, Position(1, 1, 1)));

  ReactionProfileInterpolation rpi;
  rpi.appendStructure(first, 0.0);
  ASSERT_THROW(rpi.appendStructure(second, 0.0), std::runtime_error);
}

TEST_F(ReactionProfileInterpolationTest, RPIFailure_Appending_ElementCount) {
  AtomCollection first;
  first.push_back(Atom(ElementType::He, Position(0, 0, 0)));
  first.push_back(Atom(ElementType::Xe, Position(1, 1, 1)));

  AtomCollection second;
  second.push_back(Atom(ElementType::He, Position(0, 0, 0)));
  second.push_back(Atom(ElementType::Xe, Position(1, 1, 1)));
  second.push_back(Atom(ElementType::H, Position(2, 1, 1)));

  ReactionProfileInterpolation rpi;
  rpi.appendStructure(first, 0.0);
  ASSERT_THROW(rpi.appendStructure(second, 0.0), std::runtime_error);
}

TEST_F(ReactionProfileInterpolationTest, RPIFailure_GetTS_MissingTS) {
  ReactionProfileInterpolation rpi;
  ASSERT_THROW(rpi.getCurrentTSPosition(), std::runtime_error);
}

TEST_F(ReactionProfileInterpolationTest, RPIFailure_Spline_TSisEndPoint) {
  AtomCollection first;
  first.push_back(Atom(ElementType::He, Position(0, 0, 0)));
  first.push_back(Atom(ElementType::Xe, Position(1, 1, 1)));

  AtomCollection second;
  second.push_back(Atom(ElementType::He, Position(0, 0, 0)));
  second.push_back(Atom(ElementType::Xe, Position(2, 1, 1)));

  ReactionProfileInterpolation rpi;
  rpi.appendStructure(first, 0.0);
  rpi.appendStructure(second, 0.0, true);
  ASSERT_THROW(rpi.spline(10, 3), std::runtime_error);
}

TEST_F(ReactionProfileInterpolationTest, RPIFailure_Spline_TooFewPoint) {
  AtomCollection first;
  first.push_back(Atom(ElementType::He, Position(0, 0, 0)));
  first.push_back(Atom(ElementType::Xe, Position(1, 1, 1)));

  AtomCollection second;
  second.push_back(Atom(ElementType::He, Position(0, 0, 0)));
  second.push_back(Atom(ElementType::Xe, Position(2, 1, 1)));

  AtomCollection third;
  third.push_back(Atom(ElementType::He, Position(0, 0, 0)));
  third.push_back(Atom(ElementType::Xe, Position(3, 1, 1)));

  ReactionProfileInterpolation rpi;
  rpi.appendStructure(first, 0.0);
  rpi.appendStructure(second, 0.0, true);
  rpi.appendStructure(third, 0.0);
  ASSERT_THROW(rpi.spline(2, 3), std::runtime_error);
}

TEST_F(ReactionProfileInterpolationTest, RPIFailure_Spline_TooFewPointNoTS) {
  AtomCollection first;
  first.push_back(Atom(ElementType::He, Position(0, 0, 0)));
  first.push_back(Atom(ElementType::Xe, Position(1, 1, 1)));

  AtomCollection second;
  second.push_back(Atom(ElementType::He, Position(0, 0, 0)));
  second.push_back(Atom(ElementType::Xe, Position(2, 1, 1)));

  ReactionProfileInterpolation rpi;
  rpi.appendStructure(first, 0.0);
  rpi.appendStructure(second, 0.0);
  ASSERT_THROW(rpi.spline(1, 3), std::runtime_error);
}

TEST_F(ReactionProfileInterpolationTest, SplineFailure_Eval_OutsiedInterval) {
  TrajectorySpline spline({}, {}, {});
  ASSERT_THROW(spline.evaluate(-0.001, 3), std::runtime_error);
}

TEST_F(ReactionProfileInterpolationTest, Fitting) {
  struct stat buffer;
  if (stat("Resources/dissociation.trj.xyz", &buffer) == 0) {
    ReactionProfileInterpolation rpi;
    rpi.useQuaternionFit = true;
    auto traj = MolecularTrajectoryIO::read(MolecularTrajectoryIO::format::xyz, "Resources/dissociation.trj.xyz");
    auto elements = traj.getElementTypes();
    for (int i = 0; i < traj.size(); i++) {
      rpi.appendStructure(AtomCollection(elements, traj[i]), -i);
    }
    auto spline = rpi.spline(traj.size(), 3);
    for (int i = 0; i < traj.size(); i++) {
      auto ret = spline.evaluate(double(i) / double(traj.size() - 1.0), 3);
      auto atoms = std::get<1>(ret);
      PositionCollection splinePos(atoms.getPositions());
      PositionCollection refPos(traj[i]);
      Utils::QuaternionFit fit(refPos, splinePos, elements);
      PositionCollection newPos = fit.getFittedData();
      for (unsigned int j = 0; j < refPos.rows(); j++) {
        EXPECT_NEAR((refPos.row(j) - newPos.row(j)).norm(), 0.0, 7.0e-2);
      }
    }
  }
}

TEST_F(ReactionProfileInterpolationTest, FittingAndCompression) {
  struct stat buffer;
  if (stat("Resources/dissociation.trj.xyz", &buffer) == 0) {
    ReactionProfileInterpolation rpi;
    rpi.useQuaternionFit = true;
    auto traj = MolecularTrajectoryIO::read(MolecularTrajectoryIO::format::xyz, "Resources/dissociation.trj.xyz");
    auto elements = traj.getElementTypes();
    for (int i = 0; i < traj.size(); i++) {
      rpi.appendStructure(AtomCollection(elements, traj[i]), -i);
    }
    auto spline = rpi.spline(23, 3);
    for (int i = 0; i < traj.size(); i++) {
      auto ret = spline.evaluate(double(i) / double(traj.size() - 1.0), 3);
      auto atoms = std::get<1>(ret);
      PositionCollection splinePos(atoms.getPositions());
      PositionCollection refPos(traj[i]);
      Utils::QuaternionFit fit(refPos, splinePos, elements);
      PositionCollection newPos = fit.getFittedData();
      for (unsigned int j = 0; j < refPos.rows(); j++) {
        EXPECT_NEAR((refPos.row(j) - newPos.row(j)).norm(), 0.0, 7.0e-2);
      }
    }
  }
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
