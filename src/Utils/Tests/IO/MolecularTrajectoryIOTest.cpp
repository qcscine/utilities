/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/Constants.h"
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/MolecularTrajectoryIO.h>
#include <Utils/MolecularTrajectory.h>
#include <gmock/gmock.h>
#include <boost/dll.hpp>
#include <boost/filesystem.hpp>

namespace Scine {
namespace Utils {

using namespace testing;

static void assertSameTrajectories(const MolecularTrajectory& trj, const MolecularTrajectory& trj2) {
  EXPECT_EQ(trj.size(), trj2.size());
  const auto& elementsA = trj.getElementTypes();
  const auto& elementsB = trj2.getElementTypes();
  const auto& residueInfoA = trj.getResidues();
  const auto& residueInfoB = trj.getResidues();
  for (unsigned int i = 0; i < elementsA.size(); ++i) {
    EXPECT_EQ(elementsA[i], elementsB[i]);
    if (not residueInfoA.empty()) {
      EXPECT_EQ(residueInfoA[i], residueInfoB[i]);
    }
  }
  for (int i = 0; i < trj.size(); ++i) {
    const double diff = (trj.at(i) - trj2.at(i)).array().abs().maxCoeff();
    EXPECT_NEAR(0.0, diff, 1e-4);
  }
}

TEST(AMolecularTrajectoryIOTest, CanReadMolecularTrajectory) {
  const boost::filesystem::path pathToResource = boost::dll::program_location().parent_path() / "Resources/IO_traj.xyz";
  ASSERT_NO_THROW(MolecularTrajectoryIO::read(MolecularTrajectoryIO::format::xyz, pathToResource.string()));
}

TEST(AMolecularTrajectoryIOTest, CanReadPdbTrajectory) {
  const boost::filesystem::path pathToResource =
      boost::dll::program_location().parent_path() / "Resources/trajectory.pdb";
  ASSERT_NO_THROW(MolecularTrajectoryIO::read(MolecularTrajectoryIO::format::pdb, pathToResource.string()));
  const auto trj = MolecularTrajectoryIO::read(MolecularTrajectoryIO::format::pdb, pathToResource.string());
  ElementTypeCollection referenceElements{
      ElementType::C, ElementType::C, ElementType::O, ElementType::H, ElementType::H, ElementType::H,
      ElementType::H, ElementType::H, ElementType::H, ElementType::O, ElementType::H, ElementType::H,
      ElementType::O, ElementType::H, ElementType::H, ElementType::O, ElementType::H, ElementType::H,
      ElementType::O, ElementType::H, ElementType::H, ElementType::O, ElementType::H, ElementType::H};
  const auto& trjElements = trj.getElementTypes();
  for (unsigned int i = 0; i < referenceElements.size(); ++i) {
    ASSERT_EQ(trjElements[i], referenceElements[i]);
  }
  ASSERT_EQ(trj.molecularSize(), 2520);
  ASSERT_EQ(trj.size(), 6);
  const auto& zeroPositions = trj[0];
  ASSERT_NEAR(zeroPositions(0, 0), 21.588 * Constants::bohr_per_angstrom, 1e-3);
  ASSERT_NEAR(zeroPositions(0, 1), 21.544 * Constants::bohr_per_angstrom, 1e-3);
  ASSERT_NEAR(zeroPositions(0, 2), 23.151 * Constants::bohr_per_angstrom, 1e-3);
}

TEST(AMolecularTrajectoryIOTest, ReadAndWriteTrajectory) {
  const boost::filesystem::path pathToResource =
      boost::dll::program_location().parent_path() / "Resources/trajectory.pdb";
  const auto trj = MolecularTrajectoryIO::read(MolecularTrajectoryIO::format::pdb, pathToResource.string());
  const std::string temporaryTrajectoryFile = "TemporaryTrajectoryFile.pdb";
  MolecularTrajectoryIO::write(MolecularTrajectoryIO::format::pdb, temporaryTrajectoryFile, trj);
  const auto trj2 = MolecularTrajectoryIO::read(MolecularTrajectoryIO::format::pdb, temporaryTrajectoryFile);
  assertSameTrajectories(trj, trj2);
  std::remove(temporaryTrajectoryFile.c_str());
}

TEST(AMolecularTrajectoryIOTest, ReadAndWriteTrajectoriesWithLargeCoordinates) {
  const boost::filesystem::path pathToResource =
      boost::dll::program_location().parent_path() / "Resources/large_coordinates.xyz";
  const auto trj = MolecularTrajectoryIO::read(MolecularTrajectoryIO::format::xyz, pathToResource.string());
  const std::string temporaryTrajectoryFile = "TemporaryTrajectoryFile.xyz";
  MolecularTrajectoryIO::write(MolecularTrajectoryIO::format::xyz, temporaryTrajectoryFile, trj);
  const auto trj2 = MolecularTrajectoryIO::read(MolecularTrajectoryIO::format::xyz, temporaryTrajectoryFile);
  assertSameTrajectories(trj, trj2);
  std::remove(temporaryTrajectoryFile.c_str());
}

} // namespace Utils
} // namespace Scine
