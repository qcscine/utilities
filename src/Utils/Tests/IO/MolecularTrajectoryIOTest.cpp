/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/IO/MolecularTrajectoryIO.h>
#include <Utils/MolecularTrajectory.h>
#include <gmock/gmock.h>

namespace Scine {
namespace Utils {

using namespace testing;

TEST(AMolecularTrajectoryIOTest, CanReadMolecularTrajectory) {
  struct stat buffer;
  if (stat("Resources/IO_traj.xyz", &buffer) == 0) {
    ASSERT_NO_THROW(MolecularTrajectoryIO::read(MolecularTrajectoryIO::format::xyz, "Resources/IO_traj.xyz"));
  }
}

TEST(AMolecularTrajectoryIOTest, CanReadPdbTrajectory) {
  struct stat buffer;
  if (stat("Resources/trajectory.pdb", &buffer) == 0) {
    ASSERT_NO_THROW(MolecularTrajectoryIO::read(MolecularTrajectoryIO::format::pdb, "Resources/trajectory.pdb"));
    const auto trj = MolecularTrajectoryIO::read(MolecularTrajectoryIO::format::pdb, "Resources/trajectory.pdb");
    ElementTypeCollection referenceElements{ElementType::C, ElementType::C, ElementType::O,
                                            ElementType::H, ElementType::H, ElementType::H,
                                            ElementType::H, ElementType::H, ElementType::H};
    ASSERT_EQ(trj.getElementTypes(), referenceElements);
  }
}

} // namespace Utils
} // namespace Scine
