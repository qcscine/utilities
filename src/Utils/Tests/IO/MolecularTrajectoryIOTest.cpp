/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
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
  if (stat("test_resources/traj.xyz", &buffer) == 0)
    ASSERT_NO_THROW(MolecularTrajectoryIO::read(MolecularTrajectoryIO::format::xyz, "test_resources/traj.xyz"));
}

} // namespace Utils
} // namespace Scine
