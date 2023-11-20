/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Constants.h"
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

TEST(ConstantsTest, HartreeToKCalPerMol) {
  double hartree = 1;
  double kCalPerMol = toKCalPerMol(Hartree(hartree));

  ASSERT_NEAR(kCalPerMol, 627.509, 1e-3);
}

TEST(ConstantsTest, AngstromToBohr) {
  double angstrom = 1;
  double bohr = toBohr(Angstrom(angstrom));

  ASSERT_NEAR(bohr, 1.88973, 1e-5);
}

} // namespace Tests
} /* namespace Utils */
} /* namespace Scine */
