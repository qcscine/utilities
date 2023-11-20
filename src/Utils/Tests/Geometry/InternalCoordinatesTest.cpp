/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Geometry.h"
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <gmock/gmock.h>
#include <Eigen/Geometry>

using namespace testing;
namespace Scine {
namespace Utils {
using namespace Geometry;
namespace Tests {

/**
 * @class Scine::Utils::Tests::InternalCoordinatesTest InternalCoordinatesTest.cpp
 * @brief Comprises tests for the class Scine::Utils::InternalCoordinates.
 * @test
 */
TEST(InternalCoordinatesTest, Formaldehyde) {
  auto elements = ElementTypeCollection{ElementType::O, ElementType::C, ElementType::H, ElementType::H};
  PositionCollection positions(4, 3);
  // clang-format off
  positions(0,0) = 0.0; positions(0,1) =  0.000000; positions(0,2) = -0.537500;
  positions(1,0) = 0.0; positions(1,1) =  0.000000; positions(1,2) =  0.662500;
  positions(2,0) = 0.0; positions(2,1) =  0.866025; positions(2,2) = -1.037500;
  positions(3,0) = 0.0; positions(3,1) = -0.866025; positions(3,2) = -1.037500;
  // clang-format on
  AtomCollection atoms(elements, positions);
  InternalCoordinates ic(atoms);
  auto internal = ic.coordinatesToInternal(positions);
  auto backtransformed = ic.coordinatesToCartesian(internal);
  for (unsigned int i = 0; i < 4; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      EXPECT_NEAR(positions(i, j), backtransformed(i, j), 1e-9);
    }
  }
}

TEST(InternalCoordinatesTest, Water) {
  auto elements = ElementTypeCollection{ElementType::O, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  // clang-format off
  positions(0,0) = 0.0; positions(0,1) = 0.0; positions(0,2) = 0.0;
  positions(1,0) = 0.0; positions(1,1) = 0.0; positions(1,2) = 1.0;
  positions(2,0) = 0.0; positions(2,1) = 1.0; positions(2,2) = 0.0;
  // clang-format on
  AtomCollection atoms(elements, positions);
  InternalCoordinates ic(atoms);
  ASSERT_EQ(3, ic.coordinatesToInternal(positions).size());
}

TEST(InternalCoordinatesTest, HessianGuessForSmallMolecule) {
  auto elements = ElementTypeCollection{ElementType::O, ElementType::H, ElementType::H};
  Eigen::MatrixX3d positions = Eigen::MatrixX3d::Zero(3, 3);
  // clang-format off
  positions(0,0) = 0.0; positions(0,1) = 0.0; positions(0,2) = 0.0;
  positions(1,0) = 0.0; positions(1,1) = 0.0; positions(1,2) = 1.0;
  positions(2,0) = 0.0; positions(2,1) = 1.0; positions(2,2) = 0.0;
  // clang-format on
  AtomCollection atoms(elements, positions);
  InternalCoordinates ic(atoms);
  auto guess = ic.hessianGuess();
  auto invguess = ic.inverseHessianGuess();
  Eigen::MatrixXd ref = Eigen::MatrixXd::Identity(3, 3);
  ASSERT_EQ((guess - ref).sum(), 0.0);
  ASSERT_EQ((invguess - ref).sum(), 0.0);
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
