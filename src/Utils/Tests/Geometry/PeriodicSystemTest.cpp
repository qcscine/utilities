/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Bonds/BondDetector.h"
#include "Utils/Geometry/GeometryUtilities.h"
#include "Utils/Geometry/PeriodicSystem.h"
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

/**
 * @class PeriodicSystemTest PeriodicSystemTest.cpp
 * @brief Comprises tests for the class Scine::Utils::PeriodicSystem.
 * @test
 */
class PeriodicSystemTest : public Test {
 public:
  ElementTypeCollection randomElements;
  PositionCollection randomPositions;
  PeriodicBoundaries pbc = PeriodicBoundaries(Eigen::Matrix3d::Identity());

  // PeriodicSystem shifts center of mass to center of cell -> calculate shift vector to apply back
  void shiftBack(const AtomCollection& originalAtoms, AtomCollection& shiftedAtoms) {
    auto com = Geometry::Properties::getCenterOfMass(originalAtoms);
    Position centerOfCell = Position::Constant(0.5);
    pbc.transformInPlace(centerOfCell);
    Displacement shift = com - centerOfCell;
    auto pos = shiftedAtoms.getPositions();
    Geometry::Manipulations::translatePositionsInPlace(pos, shift);
    shiftedAtoms.setPositions(pos);
  }

 private:
  void SetUp() override {
    randomElements = ElementTypeCollection{ElementType::Am, ElementType::Ce, ElementType::Ca};
    randomPositions = Eigen::MatrixX3d::Random(3, 3);
    Eigen::Matrix3d randomMatrix = Eigen::Matrix3d::Random(3, 3);
    for (int i = 0; i < 3; ++i) {
      randomMatrix(i, i) = std::fabs(randomMatrix(i, i));
    }
    pbc.setCellMatrix(randomMatrix);
  }
};

TEST_F(PeriodicSystemTest, CanBeConstructed) {
  auto ps0 = PeriodicSystem(pbc);
  ASSERT_THAT(ps0.atoms.size(), 0);
  ASSERT_THAT(ps0.pbc, pbc);
  auto ps1 = PeriodicSystem(pbc, 1);
  ASSERT_THAT(ps1.atoms.size(), 1);
  std::unordered_set<unsigned> unimportant = {1};
  auto atoms = AtomCollection(randomElements, randomPositions);
  auto ps2 = PeriodicSystem(pbc, atoms, unimportant);
  ASSERT_THAT(ps2.atoms.size(), 3);
  ASSERT_THAT(ps2.atoms.getElements(), randomElements);
}

TEST_F(PeriodicSystemTest, GetsCorrectImageAtoms) {
  /*
   *    | H H H |
   *    becomes:
   *  h | H H H | h
   *  small h represent image atoms
   */
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H, ElementType::H};
  PositionCollection positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = 0.5;
  positions(1, 0) = 1.5;
  positions(2, 0) = 2.5;
  auto h3 = AtomCollection(elements, positions);
  pbc.setCellMatrix(Eigen::Matrix3d::Identity() * 3.0);

  auto ps = PeriodicSystem(pbc, h3);
  BondOrderCollection bo = BondDetector::detectBonds(h3, pbc, true);
  ASSERT_THAT(bo, ps.constructBondOrders());
  ps.constructImageAtoms(bo);
  ps.constructImageAtoms();
  auto atoms = ps.getAtomCollectionWithImages();
  auto map = ps.getImageAtomsMap();
  ASSERT_THAT(atoms.size(), 5);
  for (const auto& element : atoms.getElements()) {
    ASSERT_THAT(element, ElementType::H);
  }
  ASSERT_DOUBLE_EQ(atoms.getPosition(3)[0], -0.5);
  ASSERT_DOUBLE_EQ(atoms.getPosition(4)[0], 3.5);
  ASSERT_NE(map.find(3), map.end());
  ASSERT_NE(map.find(4), map.end());
  ASSERT_EQ(map.at(3), 2);
  ASSERT_EQ(map.at(4), 0);
}

TEST_F(PeriodicSystemTest, GetsCorrectBondOrdersForImages) {
  /*
   *    | H-H-H |
   *      |---| <- signals that bond from right H goes to left H
   *    becomes:
   *  h-|-H-H-H-|-h
   *      |---| <- signals that bond from right H goes to left H
   *  small h represent image atoms
   */
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H, ElementType::H};
  PositionCollection positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = 0.5;
  positions(1, 0) = 1.5;
  positions(2, 0) = 2.5;
  auto h3 = AtomCollection(elements, positions);
  pbc.setCellMatrix(Eigen::Matrix3d::Identity() * 3.0);

  auto ps = PeriodicSystem(pbc, h3);
  auto data = ps.getDataForMolassemblerInterpretation();
  auto atoms = std::get<0>(data);
  const int n = 5;
  ASSERT_THAT(atoms.size(), n);
  ASSERT_DOUBLE_EQ(atoms.getPosition(3)[0], -0.5);
  ASSERT_DOUBLE_EQ(atoms.getPosition(4)[0], 3.5);
  ASSERT_THAT(ps.getAtomCollectionWithImages(), std::get<0>(data));

  BondOrderCollection bo = std::get<1>(data);
  ASSERT_THAT(bo.getSystemSize(), n);
  for (int i = 0; i < n; ++i) {
    ASSERT_THAT(bo.getOrder(i, i), 0);
  }
  ASSERT_THAT(bo.getOrder(0, 1), 1);
  ASSERT_THAT(bo.getOrder(0, 2), 1);
  ASSERT_THAT(bo.getOrder(1, 2), 1);
  ASSERT_THAT(bo.getOrder(1, 3), 0);
  ASSERT_THAT(bo.getOrder(0, 3), 0); // real atom with own image atom
  ASSERT_THAT(bo.getOrder(0, 4), 1); // real atom with image atom
  ASSERT_THAT(bo.getOrder(1, 4), 0);
  ASSERT_THAT(bo.getOrder(2, 4), 0); // real atom with own image atom
  ASSERT_THAT(bo.getOrder(2, 3), 1); // real atom with image atom
  ASSERT_THAT(ps.getBondOrderCollectionWithImages(), bo);

  ASSERT_TRUE(std::get<2>(data).empty());
  auto map = std::get<3>(data);
  ASSERT_THAT(ps.getImageAtomsMap(), map);
  ASSERT_NE(map.find(3), map.end());
  ASSERT_NE(map.find(4), map.end());
  ASSERT_EQ(map.at(3), 2);
  ASSERT_EQ(map.at(4), 0);
}

TEST_F(PeriodicSystemTest, NoImagesWithoutPeriodicBO) {
  /*
   *      H-H         H
   *    stays:
   *      H-H         H
   */
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H, ElementType::H};
  PositionCollection positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(1, 0) = 1.0;
  positions(2, 0) = 9.0;
  auto h3 = AtomCollection(elements, positions);
  pbc.setCellMatrix(Eigen::Matrix3d::Identity() * 10.0);

  auto ps = PeriodicSystem(pbc, h3);
  auto ownBonds = BondDetector::detectBonds(h3);
  auto data = ps.getDataForMolassemblerInterpretation(ownBonds);
  auto atoms = std::get<0>(data);
  const int n = 3;
  ASSERT_THAT(atoms.size(), n);

  BondOrderCollection bo = std::get<1>(data);
  ASSERT_THAT(bo.getSystemSize(), n);
  for (int i = 0; i < n; ++i) {
    ASSERT_THAT(bo.getOrder(i, i), 0);
  }
  ASSERT_THAT(bo.getOrder(0, 1), 1);
  ASSERT_THAT(bo.getOrder(0, 2), 0);
  ASSERT_THAT(bo.getOrder(1, 2), 0);

  ASSERT_TRUE(std::get<2>(data).empty());
  ASSERT_TRUE(std::get<3>(data).empty());
}

TEST_F(PeriodicSystemTest, SolidStateIndicesWork) {
  /*
   *      H-H         H
   *    stays:
   *      H-H         H
   */
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H, ElementType::H};
  PositionCollection positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(1, 0) = 1.0;
  positions(2, 0) = 9.0;
  auto h3 = AtomCollection(elements, positions);
  pbc.setCellMatrix(Eigen::Matrix3d::Identity() * 100.0);

  std::unordered_set<unsigned> unimportant = {2};
  auto ps = PeriodicSystem(pbc, h3, unimportant);
  auto ownBonds = BondDetector::detectBonds(h3);
  auto data = ps.getDataForMolassemblerInterpretation(ownBonds);
  auto atoms = std::get<0>(data);
  // PeriodicSystem shifts center of mass to center of cell -> shift atoms back
  shiftBack(h3, atoms);
  const int n = 3;
  ASSERT_THAT(atoms.size(), n);
  ASSERT_DOUBLE_EQ(atoms.getPosition(1)[0], 1.0);
  ASSERT_DOUBLE_EQ(atoms.getPosition(2)[0], 9.0);

  BondOrderCollection bo = std::get<1>(data);
  ASSERT_THAT(bo.getSystemSize(), n);
  for (int i = 0; i < n; ++i) {
    ASSERT_THAT(bo.getOrder(i, i), 0);
  }
  ASSERT_THAT(bo.getOrder(0, 1), 1);
  ASSERT_THAT(bo.getOrder(0, 2), 0);
  ASSERT_THAT(bo.getOrder(1, 2), 0);

  ASSERT_THAT(std::get<2>(data), unimportant);
  ASSERT_TRUE(std::get<3>(data).empty());
}

TEST_F(PeriodicSystemTest, SolidStateIndicesWorkPeriodic) {
  /*
   *    | H-H-H |
   *      |---| <- signals that bond from right H goes to left H
   *    becomes:
   *  h-|-H-H-H-|-h
   *      |---| <- signals that bond from right H goes to left H
   *  small h represent image atoms
   */
  auto elements = ElementTypeCollection{ElementType::H, ElementType::H, ElementType::H};
  PositionCollection positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(0, 0) = 0.5;
  positions(1, 0) = 1.5;
  positions(2, 0) = 2.5;
  auto h3 = AtomCollection(elements, positions);
  pbc.setCellMatrix(Eigen::Matrix3d::Identity() * 3.0);

  std::unordered_set<unsigned> unimportant = {0, 1, 2};
  auto ps = PeriodicSystem(pbc, h3, unimportant);
  auto data = ps.getDataForMolassemblerInterpretation();
  auto atoms = std::get<0>(data);
  // PeriodicSystem shifts center of mass to center of cell -> shift atoms back
  shiftBack(h3, atoms);
  const int n = 5;
  ASSERT_THAT(std::get<0>(data).size(), n);
  ASSERT_DOUBLE_EQ(atoms.getPosition(3)[0], -0.5);
  ASSERT_DOUBLE_EQ(atoms.getPosition(4)[0], 3.5);

  BondOrderCollection bo = std::get<1>(data);
  ASSERT_THAT(bo.getSystemSize(), n);
  for (int i = 0; i < n; ++i) {
    ASSERT_THAT(bo.getOrder(i, i), 0);
  }
  ASSERT_THAT(bo.getOrder(0, 1), 1);
  ASSERT_THAT(bo.getOrder(0, 2), 1);
  ASSERT_THAT(bo.getOrder(1, 2), 1);
  ASSERT_THAT(bo.getOrder(1, 3), 0);
  ASSERT_THAT(bo.getOrder(0, 3), 0); // real atom with own image atom
  ASSERT_THAT(bo.getOrder(0, 4), 1); // real atom with image atom
  ASSERT_THAT(bo.getOrder(1, 4), 0);
  ASSERT_THAT(bo.getOrder(2, 4), 0); // real atom with own image atom
  ASSERT_THAT(bo.getOrder(2, 3), 1); // real atom with image atom

  ASSERT_THAT(std::get<2>(data).size(), 3);
  auto map = std::get<3>(data);
  ASSERT_NE(map.find(3), map.end());
  ASSERT_NE(map.find(4), map.end());
  ASSERT_EQ(map.at(3), 2);
  ASSERT_EQ(map.at(4), 0);
}

TEST_F(PeriodicSystemTest, WorksForSurface) {
  auto elements = ElementTypeCollection{};
  const int n = 8;
  for (int i = 0; i < n; ++i) {
    elements.push_back(ElementType::H);
  }
  PositionCollection positions = Eigen::MatrixX3d::Zero(8, 3);
  // clang-format off
  positions << 0.0, 0.0, 0.0, //
               1.0, 0.0, 0.0, //            h - h - h
               0.5, 0.5, 0.0, //          h | H   H |
               0.0, 1.0, 0.0, //            H   H   h
               1.0, 1.0, 0.0, //          h | H   H |
               1.5, 0.5, 0.0, //            H - H - h
               0.5, 1.5, 0.0, //          h   h   h
               1.5, 1.5, 0.0; //
  // clang-format on
  std::unordered_set<unsigned> solid = {0, 1, 2, 3, 4, 5, 6, 7};
  pbc.setCellMatrix(Eigen::Matrix3d::Identity() * 2.0);
  auto neighbors = Geometry::Distances::countAllNearestNeighbors(positions, pbc);
  for (const auto& nNeighbors : neighbors) {
    ASSERT_THAT(nNeighbors, 4);
  }
  PeriodicSystem ps = PeriodicSystem(pbc, elements, positions, solid);
  auto bo = ps.constructBondOrders();
  bo.setToAbsoluteValues();
  for (int i = 0; i < n; ++i) {
    int counter = 0;
    for (int j = 0; j < n; ++j) {
      if (bo.getOrder(i, j) > 0.0) {
        counter++;
      }
    }
    ASSERT_THAT(counter, 4);
  }
  auto imageAtoms = ps.getImageAtoms();
  ASSERT_THAT(imageAtoms.size(), 10);
  // PeriodicSystem shifts center of mass to center of cell -> shift atoms back
  auto hSystem2d = AtomCollection(elements, positions);
  shiftBack(hSystem2d, imageAtoms);
  // clang-format off
  std::vector<std::vector<double>> expectedImagePositions = {{-0.5, -0.5, 0.0},
                                                             { 0.5, -0.5, 0.0},
                                                             { 1.5, -0.5, 0.0},
                                                             { 0.0,  2.0, 0.0},
                                                             {-0.5,  0.5, 0.0},
                                                             { 2.0,  1.0, 0.0},
                                                             {-0.5,  1.5, 0.0},
                                                             { 0.0,  2.0, 0.0},
                                                             { 1.0,  2.0, 0.0},
                                                             { 2.0,  2.0, 0.0}};
  // clang-format on
  for (const auto& expectedImage : expectedImagePositions) {
    bool found = false;
    for (int i = 0; i < imageAtoms.size(); ++i) {
      Position pos = imageAtoms.getPosition(i);
      Position expected(expectedImage.data());
      if (Geometry::Distances::distanceSquared(pos, expected) < 1e-12) {
        found = true;
      }
    }
    ASSERT_TRUE(found);
  }
}

TEST_F(PeriodicSystemTest, WorksForSurfaceShifted) {
  auto elements = ElementTypeCollection{};
  const int n = 8;
  for (int i = 0; i < n; ++i) {
    elements.push_back(ElementType::H);
  }
  PositionCollection positions = Eigen::MatrixX3d::Zero(8, 3);
  // clang-format off
  positions << 0.0, 0.0, 0.0, //
               1.0, 0.0, 0.0, //            h - h - h
               0.5, 0.5, 0.0, //          h | H   H |
               0.0, 1.0, 0.0, //            H   H   h
               1.0, 1.0, 0.0, //          h | H   H |
               1.5, 0.5, 0.0, //            H - H - h
               0.5, 1.5, 0.0, //          h   h   h
               1.5, 1.5, 0.0; //
  // clang-format on
  auto hSystem2d = AtomCollection(elements, positions);
  Displacement shift;
  shift << -1.0, 0.0, 0.0;
  Geometry::Manipulations::translatePositionsInPlace(positions, shift);
  std::unordered_set<unsigned> solid = {0, 1, 2, 3, 4, 5, 6, 7};
  pbc.setCellMatrix(Eigen::Matrix3d::Identity() * 2.0);
  auto neighbors = Geometry::Distances::countAllNearestNeighbors(positions, pbc);
  for (const auto& nNeighbors : neighbors) {
    ASSERT_THAT(nNeighbors, 4);
  }
  PeriodicSystem ps = PeriodicSystem(pbc, elements, positions, solid);
  auto bo = ps.constructBondOrders();
  bo.setToAbsoluteValues();
  for (int i = 0; i < n; ++i) {
    int counter = 0;
    for (int j = 0; j < n; ++j) {
      if (bo.getOrder(i, j) > 0.0) {
        counter++;
      }
    }
    ASSERT_THAT(counter, 4);
  }
  auto imageAtoms = ps.getImageAtoms();
  ASSERT_THAT(imageAtoms.size(), 10);
  // PeriodicSystem shifts center of mass to center of cell -> shift atoms back
  shiftBack(hSystem2d, imageAtoms);
  // clang-format off
  std::vector<std::vector<double>> expectedImagePositions = {{-0.5, -0.5, 0.0},
                                                             { 0.5, -0.5, 0.0},
                                                             { 1.5, -0.5, 0.0},
                                                             { 0.0,  2.0, 0.0},
                                                             {-0.5,  0.5, 0.0},
                                                             { 2.0,  1.0, 0.0},
                                                             {-0.5,  1.5, 0.0},
                                                             { 0.0,  2.0, 0.0},
                                                             { 1.0,  2.0, 0.0},
                                                             { 2.0,  2.0, 0.0}};
  // clang-format on
  for (const auto& expectedImage : expectedImagePositions) {
    bool found = false;
    for (int i = 0; i < imageAtoms.size(); ++i) {
      Position position = imageAtoms.getPosition(i);
      Position expected(expectedImage.data());
      if (Geometry::Distances::distanceSquared(position, expected) < 1e-12) {
        found = true;
        break;
      }
    }
    ASSERT_TRUE(found);
  }
}

TEST_F(PeriodicSystemTest, WorksForHeterogeneous3dSystem) {
  std::stringstream stream("40\n\n"
                           "Cu                 0.00000000    1.47837000   16.72589000\n"
                           "Cu                -1.28031000    3.69594000   16.72589000\n"
                           "Cu                -2.56062000    5.91350000   16.72589000\n"
                           "Cu                 2.56062000    1.47837000   16.72589000\n"
                           "Cu                 1.28031000    3.69594000   16.72589000\n"
                           "Cu                 0.00000000    5.91350000   16.72589000\n"
                           "Cu                 5.12124000    1.47837000   16.72589000\n"
                           "Cu                 3.84093000    3.69594000   16.72589000\n"
                           "Cu                 2.56062000    5.91350000   16.72589000\n"
                           "Cu                 0.00000000    0.00000000   20.90737000\n"
                           "Cu                -1.28031000    2.21756000   20.90737000\n"
                           "Cu                -2.56062000    4.43512000   20.90737000\n"
                           "Cu                 2.56062000    0.00000000   20.90737000\n"
                           "Cu                 1.28031000    2.21756000   20.90737000\n"
                           "Cu                 0.00000000    4.43512000   20.90737000\n"
                           "Cu                 5.12124000    0.00000000   20.90737000\n"
                           "Cu                 3.84093000    2.21756000   20.90737000\n"
                           "Cu                 2.56062000    4.43512000   20.90737000\n"
                           "Cu                 1.28031000    0.73919000   18.81663000\n"
                           "Cu                 0.00000000    2.95675000   18.81663000\n"
                           "Cu                -1.28031000    5.17431000   18.81663000\n"
                           "Cu                 3.84093000    0.73919000   18.81663000\n"
                           "Cu                 2.56062000    2.95675000   18.81663000\n"
                           "Cu                 1.28031000    5.17431000   18.81663000\n"
                           "Cu                 6.40155000    0.73919000   18.81663000\n"
                           "Cu                 5.12124000    2.95675000   18.81663000\n"
                           "Cu                 3.84093000    5.17431000   18.81663000\n"
                           "C                 -0.10620215    0.01981057   23.97943495\n"
                           "H                  0.18773966    0.99901306   23.66374410\n"
                           "H                  0.57373891   -0.70302889   23.57933088\n"
                           "C                 -0.08953228   -0.03591288   25.17922595\n"
                           "H                 -0.87501139    0.57308805   25.57550921\n"
                           "H                  0.85334384    0.31897779   25.53969402\n"
                           "C                 -0.29469461   -1.49463641   25.62826279\n"
                           "H                  0.56587694   -2.07332537   25.36475623\n"
                           "H                 -1.15854282   -1.89896414   25.14330366\n"
                           "C                 -0.49514096   -1.53855029   27.15453041\n"
                           "H                 -1.37424575   -0.98689652   27.41482426\n"
                           "H                 -0.60513788   -2.55483039   27.47071932\n"
                           "H                  0.35468967   -1.10443551   27.63850587\n");
  auto structure = XyzStreamHandler::read(stream);
  pbc = PeriodicBoundaries("14.516605,14.516605,71.116552,90.000000,90.000000,120.000000");
  std::unordered_set<unsigned> solidStateIndices;
  const unsigned nSurfaceAtoms = 27;
  for (unsigned i = 0; i < nSurfaceAtoms; ++i) {
    solidStateIndices.insert(i);
  }
  auto ps = Utils::PeriodicSystem(pbc, structure, solidStateIndices);
  auto data = ps.getDataForMolassemblerInterpretation();
  auto atomsWithImages = std::get<0>(data);
  auto bonds = std::get<1>(data);
  auto unimportant = std::get<2>(data);
  auto map = std::get<3>(data);
  const int nSurfaceImages = 3 * 14;
  const int nMoleculeImages = 8;
  const int nExpectedAtomsWithImages = structure.size() + nSurfaceImages + nMoleculeImages;
  ASSERT_THAT(atomsWithImages.size(), nExpectedAtomsWithImages);
  ASSERT_THAT(bonds.getSystemSize(), nExpectedAtomsWithImages);
  ASSERT_THAT(unimportant.size(), nSurfaceAtoms);
  ASSERT_THAT(map.size(), nSurfaceImages + nMoleculeImages);
}

TEST_F(PeriodicSystemTest, MultiplicationWorksSurface) {
  std::stringstream stream("27\n\n"
                           "Cu                 0.00000000    1.47837000   16.72589000\n"
                           "Cu                -1.28031000    3.69594000   16.72589000\n"
                           "Cu                -2.56062000    5.91350000   16.72589000\n"
                           "Cu                 2.56062000    1.47837000   16.72589000\n"
                           "Cu                 1.28031000    3.69594000   16.72589000\n"
                           "Cu                 0.00000000    5.91350000   16.72589000\n"
                           "Cu                 5.12124000    1.47837000   16.72589000\n"
                           "Cu                 3.84093000    3.69594000   16.72589000\n"
                           "Cu                 2.56062000    5.91350000   16.72589000\n"
                           "Cu                 0.00000000    0.00000000   20.90737000\n"
                           "Cu                -1.28031000    2.21756000   20.90737000\n"
                           "Cu                -2.56062000    4.43512000   20.90737000\n"
                           "Cu                 2.56062000    0.00000000   20.90737000\n"
                           "Cu                 1.28031000    2.21756000   20.90737000\n"
                           "Cu                 0.00000000    4.43512000   20.90737000\n"
                           "Cu                 5.12124000    0.00000000   20.90737000\n"
                           "Cu                 3.84093000    2.21756000   20.90737000\n"
                           "Cu                 2.56062000    4.43512000   20.90737000\n"
                           "Cu                 1.28031000    0.73919000   18.81663000\n"
                           "Cu                 0.00000000    2.95675000   18.81663000\n"
                           "Cu                -1.28031000    5.17431000   18.81663000\n"
                           "Cu                 3.84093000    0.73919000   18.81663000\n"
                           "Cu                 2.56062000    2.95675000   18.81663000\n"
                           "Cu                 1.28031000    5.17431000   18.81663000\n"
                           "Cu                 6.40155000    0.73919000   18.81663000\n"
                           "Cu                 5.12124000    2.95675000   18.81663000\n"
                           "Cu                 3.84093000    5.17431000   18.81663000\n");
  auto structure = XyzStreamHandler::read(stream);
  pbc = PeriodicBoundaries("14.516605,14.516605,71.116552,90.000000,90.000000,120.000000");
  const auto originalCell = pbc.getCellMatrix();
  const unsigned nSurfaceAtoms = 27;
  std::unordered_set<unsigned> solidStateIndices;
  for (unsigned i = 0; i < nSurfaceAtoms; ++i) {
    solidStateIndices.insert(i);
  }
  auto ps = PeriodicSystem(pbc, structure, solidStateIndices);
  srand(42);
  int factor = std::rand() % 10 + 1;
  auto otherPs = ps * factor;
  ASSERT_THAT(otherPs, otherPs);
  ASSERT_THAT(otherPs.atoms.size(), structure.size() * std::pow(factor, 3));
  ASSERT_THAT(otherPs.pbc.getCellMatrix(), originalCell * factor);
  ASSERT_THAT(otherPs.solidStateAtomIndices.size(), solidStateIndices.size() * std::pow(factor, 3));
  std::vector<int> randoms;
  for (int i = 0; i < 3; ++i) {
    randoms.push_back(std::rand() % 10 + 1);
  }
  Eigen::Vector3i randomScaling = Eigen::Map<Eigen::Vector3i>(randoms.data());
  auto scaledRandomPs = ps * randomScaling;
  int product = randomScaling.array().prod();
  ASSERT_THAT(scaledRandomPs.atoms.size(), structure.size() * product);
  ASSERT_THAT(scaledRandomPs.solidStateAtomIndices.size(), solidStateIndices.size() * product);
  for (const auto& index : scaledRandomPs.solidStateAtomIndices) {
    ASSERT_THAT(scaledRandomPs.atoms.getElement(index), ElementType::Cu);
  }
  ASSERT_NEAR(scaledRandomPs.pbc.getA().norm(), ps.pbc.getA().norm() * randomScaling[0], 1e-12);
  ASSERT_NEAR(scaledRandomPs.pbc.getB().norm(), ps.pbc.getB().norm() * randomScaling[1], 1e-12);
  ASSERT_NEAR(scaledRandomPs.pbc.getC().norm(), ps.pbc.getC().norm() * randomScaling[2], 1e-12);
  ASSERT_NEAR(scaledRandomPs.pbc.getAlpha(), ps.pbc.getAlpha(), 1e-12);
  ASSERT_NEAR(scaledRandomPs.pbc.getBeta(), ps.pbc.getBeta(), 1e-12);
  ASSERT_NEAR(scaledRandomPs.pbc.getGamma(), ps.pbc.getGamma(), 1e-12);
}

TEST_F(PeriodicSystemTest, MultiplicationWorksMolecular) {
  std::stringstream stream("3\n\n"
                           "H      0.7493682000    0.0000000000    0.4424329000\n"
                           "O      0.0000000000    0.0000000000   -0.1653507000\n"
                           "H     -0.7493682000    0.0000000000    0.4424329000\n");
  auto structure = XyzStreamHandler::read(stream);
  double basisLength = std::fabs(std::rand() % 10 + 1);
  pbc = PeriodicBoundaries(basisLength);
  auto ps = PeriodicSystem(pbc, structure);
  srand(42);
  int factor = std::rand() % 10 + 1;
  auto otherPs = ps * factor;
  ps *= factor;
  ASSERT_THAT(ps, otherPs);
  ASSERT_THAT(ps.atoms.size(), structure.size() * std::pow(factor, 3));
  ASSERT_THAT(ps.pbc.getCellMatrix(), Eigen::Matrix3d::Identity() * factor * basisLength);
  ASSERT_TRUE(ps.solidStateAtomIndices.empty());
  auto rm = Eigen::Matrix3d::Random().cwiseAbs();
  Eigen::Matrix3d randomMatrix = rm;
  randomMatrix += Eigen::Matrix3d::Ones();
  std::unordered_set<unsigned> testIndices = {0, 2};
  auto randomPs = PeriodicSystem(PeriodicBoundaries(randomMatrix), structure, testIndices);
  Eigen::Vector3i scaling;
  scaling << 3, 1, 2;
  PeriodicSystem scaledRandomPs = randomPs * scaling;
  ASSERT_THAT(scaledRandomPs.atoms.size(), structure.size() * 3 * 2);
  ASSERT_THAT(scaledRandomPs.solidStateAtomIndices.size(), testIndices.size() * 3 * 2);
  ASSERT_NEAR(scaledRandomPs.pbc.getA().norm(), randomPs.pbc.getA().norm() * 3, 1e-12);
  ASSERT_NEAR(scaledRandomPs.pbc.getB().norm(), randomPs.pbc.getB().norm() * 1, 1e-12);
  ASSERT_NEAR(scaledRandomPs.pbc.getC().norm(), randomPs.pbc.getC().norm() * 2, 1e-12);
  ASSERT_NEAR(scaledRandomPs.pbc.getAlpha(), randomPs.pbc.getAlpha(), 1e-12);
  ASSERT_NEAR(scaledRandomPs.pbc.getBeta(), randomPs.pbc.getBeta(), 1e-12);
  ASSERT_NEAR(scaledRandomPs.pbc.getGamma(), randomPs.pbc.getGamma(), 1e-12);
}

TEST_F(PeriodicSystemTest, ConstructingNegativeBondOrdersWorks) {
  std::stringstream stream("40\n\n"
                           "Cu                 0.00000000    1.47837000   16.72589000\n"
                           "Cu                -1.28031000    3.69594000   16.72589000\n"
                           "Cu                -2.56062000    5.91350000   16.72589000\n"
                           "Cu                 2.56062000    1.47837000   16.72589000\n"
                           "Cu                 1.28031000    3.69594000   16.72589000\n"
                           "Cu                 0.00000000    5.91350000   16.72589000\n"
                           "Cu                 5.12124000    1.47837000   16.72589000\n"
                           "Cu                 3.84093000    3.69594000   16.72589000\n"
                           "Cu                 2.56062000    5.91350000   16.72589000\n"
                           "Cu                 0.00000000    0.00000000   20.90737000\n"
                           "Cu                -1.28031000    2.21756000   20.90737000\n"
                           "Cu                -2.56062000    4.43512000   20.90737000\n"
                           "Cu                 2.56062000    0.00000000   20.90737000\n"
                           "Cu                 1.28031000    2.21756000   20.90737000\n"
                           "Cu                 0.00000000    4.43512000   20.90737000\n"
                           "Cu                 5.12124000    0.00000000   20.90737000\n"
                           "Cu                 3.84093000    2.21756000   20.90737000\n"
                           "Cu                 2.56062000    4.43512000   20.90737000\n"
                           "Cu                 1.28031000    0.73919000   18.81663000\n"
                           "Cu                 0.00000000    2.95675000   18.81663000\n"
                           "Cu                -1.28031000    5.17431000   18.81663000\n"
                           "Cu                 3.84093000    0.73919000   18.81663000\n"
                           "Cu                 2.56062000    2.95675000   18.81663000\n"
                           "Cu                 1.28031000    5.17431000   18.81663000\n"
                           "Cu                 6.40155000    0.73919000   18.81663000\n"
                           "Cu                 5.12124000    2.95675000   18.81663000\n"
                           "Cu                 3.84093000    5.17431000   18.81663000\n"
                           "C                 -0.10620215    0.01981057   23.97943495\n"
                           "H                  0.18773966    0.99901306   23.66374410\n"
                           "H                  0.57373891   -0.70302889   23.57933088\n"
                           "C                 -0.08953228   -0.03591288   25.17922595\n"
                           "H                 -0.87501139    0.57308805   25.57550921\n"
                           "H                  0.85334384    0.31897779   25.53969402\n"
                           "C                 -0.29469461   -1.49463641   25.62826279\n"
                           "H                  0.56587694   -2.07332537   25.36475623\n"
                           "H                 -1.15854282   -1.89896414   25.14330366\n"
                           "C                 -0.49514096   -1.53855029   27.15453041\n"
                           "H                 -1.37424575   -0.98689652   27.41482426\n"
                           "H                 -0.60513788   -2.55483039   27.47071932\n"
                           "H                  0.35468967   -1.10443551   27.63850587\n");
  auto structure = XyzStreamHandler::read(stream);
  pbc = PeriodicBoundaries("14.516605,14.516605,71.116552,90.000000,90.000000,120.000000");
  std::unordered_set<unsigned> solidStateIndices;
  const unsigned nSurfaceAtoms = 27;
  for (unsigned i = 0; i < nSurfaceAtoms; ++i) {
    solidStateIndices.insert(i);
  }
  auto ps = Utils::PeriodicSystem(pbc, structure, solidStateIndices);
  auto bo = ps.constructBondOrders();
  const BondOrderCollection refBo = bo;
  bo.setToAbsoluteValues();
  ps.makeBondOrdersAcrossBoundariesNegative(bo);
  ASSERT_TRUE(refBo == bo);
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
