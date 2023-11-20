/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Math/MachineLearning/PrincipalComponentAnalysis.h>
#include <gmock/gmock.h>

namespace Scine {
namespace Utils {
using namespace MachineLearning;
using namespace testing;
namespace Tests {

/**
 * @class APrincipalComponentAnalysisTest @file PrincipalComponentAnalysisTest.cpp
 * @test
 *
 * The PCA is tested against the implementation of scikit-learn:
 * https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html (visited Jan 13, 2020)
 *
 */
class APrincipalComponentAnalysisTest : public Test {
 public:
  Eigen::MatrixXd data;
  Eigen::VectorXd correctFirstComponent;
  Eigen::VectorXd correctSecondComponent;
  Eigen::VectorXd correctThirdComponent;

 private:
  void SetUp() final {
    data.resize(7, 3);
    data(0, 0) = 0.95291535;
    data(0, 1) = 0.98657902;
    data(0, 2) = 0.63002747;
    data(1, 0) = 0.38204747;
    data(1, 1) = 0.26429569;
    data(1, 2) = 0.08035464;
    data(2, 0) = 0.98627932;
    data(2, 1) = 0.08674566;
    data(2, 2) = 0.67170195;
    data(3, 0) = 0.05343161;
    data(3, 1) = 0.25514023;
    data(3, 2) = 0.36168503;
    data(4, 0) = 0.43645506;
    data(4, 1) = 0.73115828;
    data(4, 2) = 0.06223424;
    data(5, 0) = 0.80503707;
    data(5, 1) = 0.13232257;
    data(5, 2) = 0.85751018;
    data(6, 0) = 0.98791868;
    data(6, 1) = 0.09372797;
    data(6, 2) = 0.03823866;

    correctFirstComponent.resize(3);
    correctFirstComponent(0) = -0.77620393;
    correctFirstComponent(1) = 0.10553319;
    correctFirstComponent(2) = -0.62158684;
    correctSecondComponent.resize(3);
    correctSecondComponent(0) = 0.10460237;
    correctSecondComponent(1) = 0.99378393;
    correctSecondComponent(2) = 0.03810315;
    correctThirdComponent.resize(3);
    correctThirdComponent(0) = 0.62174416;
    correctThirdComponent(1) = -0.03544364;
    correctThirdComponent(2) = -0.78241801;
  }
};

TEST_F(APrincipalComponentAnalysisTest, ComponentsAreCorrectlyObtained) {
  PrincipalComponentAnalysis pca1(data);
  auto result = pca1.calculate(3);

  auto firstComponent = result.first.col(0);
  auto secondComponent = result.first.col(1);
  auto thirdComponent = result.first.col(2);

  EXPECT_THAT(result.first.size(), Eq(9));
  EXPECT_THAT(result.second.size(), Eq(3));

  // Check that the norms are correct
  EXPECT_THAT(firstComponent.squaredNorm(), DoubleNear(correctFirstComponent.squaredNorm(), 1e-8));
  EXPECT_THAT(secondComponent.squaredNorm(), DoubleNear(correctSecondComponent.squaredNorm(), 1e-8));
  EXPECT_THAT(thirdComponent.squaredNorm(), DoubleNear(correctThirdComponent.squaredNorm(), 1e-8));

  // Check that they are also parallel
  EXPECT_THAT(std::abs(firstComponent.dot(correctFirstComponent)),
              DoubleNear(firstComponent.norm() * correctFirstComponent.norm(), 1e-6));
  EXPECT_THAT(std::abs(secondComponent.dot(correctSecondComponent)),
              DoubleNear(secondComponent.norm() * correctSecondComponent.norm(), 1e-6));
  EXPECT_THAT(std::abs(thirdComponent.dot(correctThirdComponent)),
              DoubleNear(thirdComponent.norm() * correctThirdComponent.norm(), 1e-6));

  EXPECT_THAT(result.second(2), DoubleNear(0.20055418, 1e-6));
  EXPECT_THAT(result.second(1), DoubleNear(0.33143852, 1e-6));
  EXPECT_THAT(result.second(0), DoubleNear(0.4680073, 1e-6));

  PrincipalComponentAnalysis pca2(data);
  auto resultWithOnlyOneComponent = pca2.calculate(1);

  EXPECT_THAT(resultWithOnlyOneComponent.first.size(), Eq(3));
  EXPECT_THAT(resultWithOnlyOneComponent.second.size(), Eq(1));

  EXPECT_THAT(resultWithOnlyOneComponent.first.col(0).squaredNorm(), DoubleNear(correctFirstComponent.squaredNorm(), 1e-8));
  EXPECT_THAT(std::abs(resultWithOnlyOneComponent.first.col(0).dot(correctFirstComponent)),
              DoubleNear(resultWithOnlyOneComponent.first.col(0).norm() * correctFirstComponent.norm(), 1e-6));
  EXPECT_THAT(resultWithOnlyOneComponent.second(0), DoubleNear(0.4680073, 1e-6));
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
