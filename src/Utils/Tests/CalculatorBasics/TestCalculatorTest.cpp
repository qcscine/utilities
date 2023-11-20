/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
/* Tested File */
#include "Utils/CalculatorBasics/TestCalculator.h"
/* Internal */
#include "Utils/CalculatorBasics.h"
#include "Utils/CalculatorBasics/PropertyList.h"
#include "Utils/Geometry/AtomCollection.h"
/* External */
#include <gmock/gmock.h>

using namespace testing;

namespace Scine {
namespace Utils {
namespace Tests {

TEST(TestCalculatorTests, EnergyAndGradients) {
  auto elements = ElementTypeCollection{ElementType::Cl, ElementType::C, ElementType::H,
                                        ElementType::H,  ElementType::H, ElementType::Br};
  PositionCollection startPositions = PositionCollection::Zero(6, 3);
  // clang-format off
  startPositions <<  3.7376961460e+00,  2.1020866350e-04,  4.5337439168e-02,
                    -3.8767703481e+00, -2.4803422157e-05, -1.2049608882e-01,
                    -2.3620148614e+00,  1.3238308540e+00,  1.0376490681e-01,
                    -2.3809041075e+00, -8.2773666259e-01,  9.6331578315e-01,
                    -2.3309449521e+00, -4.9652606314e-01, -1.3293307598e+00,
                    -7.4798903722e+00,  2.6536371103e-04, -1.9897114399e-01;
  // clang-format on
  AtomCollection atoms(elements, startPositions);

  TestCalculator testCalculator;
  testCalculator.setStructure(atoms);
  auto results = testCalculator.calculate();
  double energy = results.get<Utils::Property::Energy>();
  auto gradients = results.get<Property::Gradients>();

  double refEnergy = -2.91286636256899e-01;
  GradientCollection refGradients = GradientCollection::Zero(6, 3);
  // clang-format off
  refGradients << 1.8200000000e-11,   1.2419900000e-07,  -1.1075000000e-09,
                  9.0500000000e-11,  -2.0064900000e-08,   1.4420000000e-10,
                  5.0392400000e-08,  -7.8658100000e-08,   2.1288000000e-09,
                 -3.1917600000e-08,  -7.9348100000e-08,  -7.3340000000e-10,
                 -1.8650900000e-08,  -7.9893700000e-08,   7.7700000000e-10,
                  6.7200000000e-11,   1.3376590000e-07,  -1.2091000000e-09;
  // clang-format on

  EXPECT_NEAR(refEnergy, energy, 1e-13);
  for (unsigned int i = 0; i < refGradients.array().size(); i++) {
    EXPECT_NEAR(refGradients.data()[i], gradients.data()[i], 1.0e-12);
  }
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
