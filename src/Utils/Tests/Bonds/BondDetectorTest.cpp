/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Bonds/BondDetector.h"
#include "Utils/Bonds/BondOrderCollection.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/IO/ChemicalFileFormats/MolStreamHandler.h"
#include <gmock/gmock.h>
#include <Eigen/SparseCore>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

TEST(BondDetectorTest, Empty) {
  auto emptyAc = AtomCollection();
  BondOrderCollection emptyBO = Utils::BondDetector::detectBonds(emptyAc);
  ASSERT_EQ(emptyBO.getSystemSize(), 0);
  ASSERT_TRUE(emptyBO.empty());

  int size = 3;
  emptyAc = AtomCollection(size);
  emptyBO = Utils::BondDetector::detectBonds(emptyAc);
  ASSERT_EQ(emptyBO.getSystemSize(), size);
  ASSERT_EQ(emptyBO.getOrder(0, 0), 0);
  ASSERT_EQ(emptyBO.getMatrix().diagonal(), Eigen::VectorXd::Zero(size));
  // ac init with int sets all positions to 0.0 -> all bonds
  for (int i = 0; i < size - 1; ++i) {
    for (int j = i + 1; j < size; ++j) {
      ASSERT_EQ(emptyBO.getOrder(i, j), 1);
    }
  }
}

TEST(BondDetectorTest, H2Correct) {
  PositionCollection positions = Eigen::MatrixX3d::Zero(2, 3);
  positions(1, 0) = 0.6;
  Utils::ElementTypeCollection elements;
  elements.push_back(Utils::ElementType::H);
  elements.push_back(Utils::ElementType::H);

  auto h2 = AtomCollection(elements, positions);
  BondOrderCollection bo = Utils::BondDetector::detectBonds(h2);
  ASSERT_EQ(bo.getSystemSize(), 2);
  ASSERT_EQ(bo.getOrder(0, 0), 0);
  ASSERT_EQ(bo.getOrder(0, 1), 1);
  ASSERT_EQ(bo.getOrder(1, 0), 1);
  ASSERT_EQ(bo.getOrder(1, 1), 0);
}

TEST(BondDetectorTest, H2CorrectPeriodic) {
  PositionCollection positions = Eigen::MatrixX3d::Zero(2, 3);
  positions(1, 0) = 0.6;
  Utils::ElementTypeCollection elements;
  elements.push_back(Utils::ElementType::H);
  elements.push_back(Utils::ElementType::H);

  auto h2 = AtomCollection(elements, positions);
  PeriodicBoundaries pbc(Eigen::Matrix3d::Identity() * 10.0);
  BondOrderCollection bo = Utils::BondDetector::detectBonds(h2, pbc);
  ASSERT_EQ(bo.getSystemSize(), 2);
  ASSERT_EQ(bo.getOrder(0, 0), 0);
  ASSERT_EQ(bo.getOrder(0, 1), 1);
  ASSERT_EQ(bo.getOrder(1, 0), 1);
  ASSERT_EQ(bo.getOrder(1, 1), 0);
}

TEST(BondDetectorTest, CorrectNegativePeriodicOrders) {
  PositionCollection positions = Eigen::MatrixX3d::Zero(3, 3);
  positions(1, 0) = 1.0;
  positions(2, 0) = 9.0;
  Utils::ElementTypeCollection elements;
  elements.push_back(Utils::ElementType::H);
  elements.push_back(Utils::ElementType::H);
  elements.push_back(Utils::ElementType::H);

  auto h3 = AtomCollection(elements, positions);
  PeriodicBoundaries pbc(Eigen::Matrix3d::Identity() * 10.0);
  BondOrderCollection bo = Utils::BondDetector::detectBonds(h3, pbc, true);
  ASSERT_EQ(bo.getSystemSize(), 3);
  ASSERT_EQ(bo.getOrder(0, 0), 0);
  ASSERT_EQ(bo.getOrder(0, 1), 1);
  ASSERT_EQ(bo.getOrder(0, 2), -1);
  ASSERT_EQ(bo.getOrder(2, 0), -1);
  ASSERT_EQ(bo.getOrder(1, 0), 1);
  ASSERT_EQ(bo.getOrder(1, 1), 0);
  ASSERT_EQ(bo.getOrder(1, 2), 0);
  ASSERT_EQ(bo.getOrder(2, 1), 0);
  ASSERT_EQ(bo.getOrder(2, 2), 0);
}

TEST(BondDetectorTest, H2Edge) {
  PositionCollection positions = Eigen::MatrixX3d::Zero(2, 3);
  // limit with 0.4 Angstrom tolerance for H2 is 1.625164 Bohr
  positions(1, 0) = 1.62516;
  Utils::ElementTypeCollection elements;
  elements.push_back(Utils::ElementType::H);
  elements.push_back(Utils::ElementType::H);

  auto h2 = AtomCollection(elements, positions);
  BondOrderCollection bo = Utils::BondDetector::detectBonds(h2);
  ASSERT_EQ(bo.getSystemSize(), 2);
  ASSERT_EQ(bo.getOrder(0, 0), 0);
  ASSERT_EQ(bo.getOrder(0, 1), 1);
  ASSERT_EQ(bo.getOrder(1, 0), 1);
  ASSERT_EQ(bo.getOrder(1, 1), 0);
}

TEST(BondDetectorTest, H2Wrong) {
  PositionCollection positions = Eigen::MatrixX3d::Zero(2, 3);
  positions(1, 0) = 1.8;
  Utils::ElementTypeCollection elements;
  elements.push_back(Utils::ElementType::H);
  elements.push_back(Utils::ElementType::H);

  auto h2 = AtomCollection(elements, positions);
  BondOrderCollection bo = Utils::BondDetector::detectBonds(h2);
  ASSERT_EQ(bo.getSystemSize(), 2);
  ASSERT_TRUE(bo.empty());

  // but bond with VdW
  bo = Utils::BondDetector::detectBonds(h2, true);
  ASSERT_EQ(bo.getSystemSize(), 2);
  ASSERT_FALSE(bo.empty());
  ASSERT_EQ(bo.getOrder(0, 0), 0);
  ASSERT_EQ(bo.getOrder(0, 1), 1);
  ASSERT_EQ(bo.getOrder(1, 0), 1);
  ASSERT_EQ(bo.getOrder(1, 1), 0);
}

TEST(BondDetectorTest, HCNCorrect) {
  std::string hcn = "hcn\n"
                    "opt with PBE def2-svp then PyMOL conversion\n"
                    "\n"
                    "  3  2  0  0  0  0  0  0  0  0999 V2000\n"
                    "    0.0000    0.0000    1.0432 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                    "   -0.0000   -0.0000    2.2104 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
                    "   -0.0000   -0.0000   -0.0440 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                    "  1  3  1  0  0  0\n"
                    "  1  2  3  0  0  0\n"
                    "M  END";
  std::stringstream in(hcn);
  auto data = MolStreamHandler::read(in);
  AtomCollection atoms = data.first;
  BondOrderCollection refBO = data.second;
  BondOrderCollection bo = Utils::BondDetector::detectBonds(atoms);
  ASSERT_EQ(bo.getSystemSize(), refBO.getSystemSize());
  for (int i = 0; i < bo.getSystemSize(); ++i) {
    for (int j = 0; j < bo.getSystemSize(); ++j) {
      if (refBO.getOrder(i, j) > 0.1)
        ASSERT_EQ(bo.getOrder(i, j), 1);
      else
        ASSERT_EQ(bo.getOrder(i, j), 0);
    }
  }
}

TEST(BondDetectorTest, FerroceneCorrect) {
  std::string ferrocene = "ferrocene\n"
                          "opt with PBE def2-svp then PyMOL conversion\n"
                          "\n"
                          " 21 30  0  0  0  0  0  0  0  0999 V2000\n"
                          "    1.6683    0.1221   -0.2724 Fe  0  0  0  0  0  0  0  0  0  0  0  0\n"
                          "    2.7959   -0.8838   -1.6423 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                          "    2.9122    0.5382    1.2889 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                          "    2.2704    0.3315   -2.2088 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                          "    1.6843   -1.6902   -1.2096 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                          "    1.8014   -0.2675    1.7247 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                          "    2.3853    1.7521    0.7207 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                          "    0.8338    0.2759   -2.1264 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                          "    0.4716   -0.9737   -1.5086 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                          "    0.5878    0.4480    1.4266 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                          "    0.9490    1.6962    0.8057 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                          "    3.8579   -1.1440   -1.5472 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                          "    3.9736    0.2704    1.3665 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                          "    2.8612    1.1592   -2.6217 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                          "    1.7502   -2.6742   -0.7290 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                          "    1.8680   -1.2573    2.1941 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                          "    2.9744    2.5725    0.2912 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                          "    0.1375    1.0537   -2.4644 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                          "   -0.5490   -1.3141   -1.2928 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                          "   -0.4329    0.1007    1.6305 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                          "    0.2515    2.4670    0.4532 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                          "  1  2  1  0  0  0  0\n"
                          "  1  3  1  0  0  0  0\n"
                          "  1  4  1  0  0  0  0\n"
                          "  1  5  1  0  0  0  0\n"
                          "  1  6  1  0  0  0  0\n"
                          "  1  7  1  0  0  0  0\n"
                          "  1  8  1  0  0  0  0\n"
                          "  1  9  1  0  0  0  0\n"
                          "  1 10  1  0  0  0  0\n"
                          "  1 11  1  0  0  0  0\n"
                          "  2  4  4  0  0  0  0\n"
                          "  2  5  4  0  0  0  0\n"
                          "  2 12  1  0  0  0  0\n"
                          "  3  6  4  0  0  0  0\n"
                          "  3  7  4  0  0  0  0\n"
                          "  3 13  1  0  0  0  0\n"
                          "  4  8  4  0  0  0  0\n"
                          "  4 14  1  0  0  0  0\n"
                          "  5  9  4  0  0  0  0\n"
                          "  5 15  1  0  0  0  0\n"
                          "  6 10  4  0  0  0  0\n"
                          "  6 16  1  0  0  0  0\n"
                          "  7 11  4  0  0  0  0\n"
                          "  7 17  1  0  0  0  0\n"
                          "  8  9  4  0  0  0  0\n"
                          "  8 18  1  0  0  0  0\n"
                          "  9 19  1  0  0  0  0\n"
                          " 10 11  4  0  0  0  0\n"
                          " 10 20  1  0  0  0  0\n"
                          " 11 21  1  0  0  0  0\n"
                          "M  END\n";
  std::stringstream in(ferrocene);
  auto data = MolStreamHandler::read(in);
  AtomCollection atoms = data.first;
  BondOrderCollection refBO = data.second;
  BondOrderCollection bo = Utils::BondDetector::detectBonds(atoms);
  ASSERT_EQ(bo.getSystemSize(), refBO.getSystemSize());
  for (int i = 0; i < bo.getSystemSize(); ++i) {
    for (int j = 0; j < bo.getSystemSize(); ++j) {
      if (refBO.getOrder(i, j) > 0.1)
        ASSERT_EQ(bo.getOrder(i, j), 1);
      else
        ASSERT_EQ(bo.getOrder(i, j), 0);
    }
  }
}

TEST(BondDetectorTest, ClafenCorrect) {
  std::string clafen = "clafen\n"
                       "opt with PBE def2-svp then PyMOL conversion\n"
                       "\n"
                       " 29 29  0  0  0  0  0  0  0  0999 V2000\n"
                       "   -1.6254   -0.4809    0.9857 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "   -2.0566    1.0662   -1.0163 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "   -3.0289   -0.4027    1.3133 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "   -3.5992    0.9552    0.9105 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "   -3.4630    1.2028   -0.5940 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "    1.3718   -0.8498   -1.2022 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "    1.6782   -2.0848   -0.3506 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "    2.9217   -3.1107   -1.1608 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "    0.4035    0.0619   -0.5831 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "    0.8694    1.2280    0.1635 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "    1.2969    0.9118    1.6012 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "    1.7625    2.4306    2.4590 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "   -1.2236   -0.3683   -0.6156 P   0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "   -1.4675   -1.6017   -1.4395 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "   -1.9179    1.3057   -2.0049 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "   -3.5744   -1.2323    0.8085 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "   -3.0980   -0.5654    2.4072 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "   -4.6699    1.0062    1.1982 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "   -3.0647    1.7510    1.4708 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "   -3.7982    2.2327   -0.8362 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "   -4.1361    0.5047   -1.1483 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "    2.3056   -0.2819   -1.3954 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "    0.9813   -1.1960   -2.1799 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "    2.0722   -1.8173    0.6485 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "    0.7642   -2.6983   -0.2335 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "    1.7167    1.6976   -0.3817 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "    0.0491    1.9742    0.1901 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "    0.4613    0.4461    2.1570 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "    2.1781    0.2410    1.6387 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                       "  1  3  1  0  0  0  0\n"
                       "  1 13  1  0  0  0  0\n"
                       "  2  5  1  0  0  0  0\n"
                       "  2 13  1  0  0  0  0\n"
                       "  2 15  1  0  0  0  0\n"
                       "  3  4  1  0  0  0  0\n"
                       "  3 16  1  0  0  0  0\n"
                       "  3 17  1  0  0  0  0\n"
                       "  4  5  1  0  0  0  0\n"
                       "  4 18  1  0  0  0  0\n"
                       "  4 19  1  0  0  0  0\n"
                       "  5 20  1  0  0  0  0\n"
                       "  5 21  1  0  0  0  0\n"
                       "  6  7  1  0  0  0  0\n"
                       "  6  9  1  0  0  0  0\n"
                       "  6 22  1  0  0  0  0\n"
                       "  6 23  1  0  0  0  0\n"
                       "  7  8  1  0  0  0  0\n"
                       "  7 24  1  0  0  0  0\n"
                       "  7 25  1  0  0  0  0\n"
                       "  9 10  1  0  0  0  0\n"
                       "  9 13  1  0  0  0  0\n"
                       " 10 11  1  0  0  0  0\n"
                       " 10 26  1  0  0  0  0\n"
                       " 10 27  1  0  0  0  0\n"
                       " 11 12  1  0  0  0  0\n"
                       " 11 28  1  0  0  0  0\n"
                       " 11 29  1  0  0  0  0\n"
                       " 13 14  2  0  0  0  0\n"
                       "M  END";
  std::stringstream in(clafen);
  auto data = MolStreamHandler::read(in);
  AtomCollection atoms = data.first;
  BondOrderCollection refBO = data.second;
  BondOrderCollection bo = Utils::BondDetector::detectBonds(atoms);
  ASSERT_EQ(bo.getSystemSize(), refBO.getSystemSize());
  for (int i = 0; i < bo.getSystemSize(); ++i) {
    for (int j = 0; j < bo.getSystemSize(); ++j) {
      if (refBO.getOrder(i, j) > 0.1)
        ASSERT_EQ(bo.getOrder(i, j), 1);
      else
        ASSERT_EQ(bo.getOrder(i, j), 0);
    }
  }
}

TEST(BondDetectorTest, DimethylButaneCorrect) {
  std::string dimethylbutane = "C6H14 Butane, 2,2-dimethyl- 75832\n"
                               "##CCCBDB10101603:55\n"
                               "Geometry Optimized at B3LYP/TZVP\n"
                               " 20 19  0  0  0  0  0  0  0  0999 V2000\n"
                               "    1.9046   -1.0745    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                               "    1.1702    0.2696    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                               "   -0.3765    0.2200    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                               "   -0.8983    1.6682    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                               "   -0.8983   -0.4979    1.2565 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                               "   -0.8983   -0.4979   -1.2565 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                               "   -1.9910    1.6914    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                               "    2.9849   -0.9134    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                               "   -0.5247   -0.0194    2.1656 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                               "   -0.5247   -0.0194   -2.1656 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                               "   -0.5966   -1.5471    1.2808 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                               "   -0.5966   -1.5471   -1.2808 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                               "   -1.9906   -0.4711    1.2894 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                               "   -1.9906   -0.4711   -1.2894 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                               "    1.6651   -1.6718    0.8821 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                               "    1.6651   -1.6718   -0.8821 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                               "    1.4911    0.8456    0.8748 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                               "    1.4911    0.8456   -0.8748 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                               "   -0.5520    2.2123    0.8829 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                               "   -0.5520    2.2123   -0.8829 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                               "  1  2  1  0  0  0\n"
                               "  1  8  1  0  0  0\n"
                               "  1 15  1  0  0  0\n"
                               "  1 16  1  0  0  0\n"
                               "  2  3  1  0  0  0\n"
                               "  2 17  1  0  0  0\n"
                               "  2 18  1  0  0  0\n"
                               "  3  4  1  0  0  0\n"
                               "  3  5  1  0  0  0\n"
                               "  3  6  1  0  0  0\n"
                               "  4  7  1  0  0  0\n"
                               "  4 19  1  0  0  0\n"
                               "  4 20  1  0  0  0\n"
                               "  5  9  1  0  0  0\n"
                               "  5 11  1  0  0  0\n"
                               "  5 13  1  0  0  0\n"
                               "  6 10  1  0  0  0\n"
                               "  6 12  1  0  0  0\n"
                               "  6 14  1  0  0  0\n"
                               "M  END";
  std::stringstream in(dimethylbutane);
  auto data = MolStreamHandler::read(in);
  AtomCollection atoms = data.first;
  BondOrderCollection refBO = data.second;
  BondOrderCollection bo = Utils::BondDetector::detectBonds(atoms);
  ASSERT_EQ(bo.getSystemSize(), refBO.getSystemSize());
  for (int i = 0; i < bo.getSystemSize(); ++i) {
    for (int j = 0; j < bo.getSystemSize(); ++j) {
      if (refBO.getOrder(i, j) > 0.1)
        ASSERT_EQ(bo.getOrder(i, j), 1);
      else
        ASSERT_EQ(bo.getOrder(i, j), 0);
    }
  }
}

TEST(BondDetectorTest, PeriodicCuCorrect) {
  std::string cu = "Cu\n"
                   "Materialsproject, conversion with PyMOL\n"
                   "\n"
                   "  4  6  0  0  0  0  0  0  0  0999 V2000\n"
                   "    0.0000    0.0000    0.0000 Cu  0  0  0  0  0  0  0  0  0  0  0  0\n"
                   "   -0.0000    1.8106    1.8106 Cu  0  0  0  0  0  0  0  0  0  0  0  0\n"
                   "    1.8106   -0.0000    1.8106 Cu  0  0  0  0  0  0  0  0  0  0  0  0\n"
                   "    1.8106    1.8106    0.0000 Cu  0  0  0  0  0  0  0  0  0  0  0  0\n"
                   "  2  4  1  0  0  0  0\n"
                   "  3  4  1  0  0  0  0\n"
                   "  1  4  1  0  0  0  0\n"
                   "  2  3  1  0  0  0  0\n"
                   "  1  2  1  0  0  0  0\n"
                   "  1  3  1  0  0  0  0\n"
                   "M  END";
  std::stringstream in(cu);
  auto data = MolStreamHandler::read(in);
  AtomCollection atoms = data.first;
  BondOrderCollection refBO = data.second;
  Eigen::Vector3d lengths;
  Eigen::Vector3d angles;
  lengths << 2.561, 2.561, 2.561;
  angles << 90.0, 90.0, 90.0;
  auto pbc = PeriodicBoundaries(lengths, angles, false);
  BondOrderCollection bo = Utils::BondDetector::detectBonds(atoms, pbc);
  ASSERT_EQ(bo.getSystemSize(), refBO.getSystemSize());
  ASSERT_EQ(bo.getMatrix().diagonal(), Eigen::VectorXd::Zero(bo.getSystemSize()));
  for (int i = 0; i < bo.getSystemSize() - 1; ++i) {
    for (int j = i + 1; j < bo.getSystemSize(); ++j) {
      ASSERT_EQ(bo.getOrder(i, j), 1);
    }
  }
}

TEST(BondDetectorTest, PeriodicCoOCorrect) {
  std::string coo = "coo\n"
                    "MaterialsProject and Conversion with PyMOL after setting bonds\n"
                    "\n"
                    "  4  6  0  0  0  0  0  0  0  0999 V2000\n"
                    "    1.6485    0.9518    4.6167 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
                    "   -0.0000    1.9036    1.9896 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
                    "    1.6485    0.9518    2.6078 Co  0  0  0  0  0  0  0  0  0  0  0  0\n"
                    "   -0.0000    1.9036    5.2349 Co  0  0  0  0  0  0  0  0  0  0  0  0\n"
                    "  1  3  1  0  0  0  0\n"
                    "  1  4  1  0  0  0  0\n"
                    "  2  3  1  0  0  0  0\n"
                    "  1  4  1  0  0  0  0\n"
                    "  1  4  1  0  0  0  0\n"
                    "  2  4  1  0  0  0  0\n"
                    "M  END";
  std::stringstream in(coo);
  auto data = MolStreamHandler::read(in);
  AtomCollection atoms = data.first;
  BondOrderCollection refBO = data.second;
  Eigen::Vector3d lengths;
  Eigen::Vector3d angles;
  double a = 3.297078325;
  double b = 3.29707832413;
  double c = 5.254212582;
  double alpha = 90.0;
  double beta = 90.0;
  double gamma = 120.0;
  lengths << a, b, c;
  angles << alpha, beta, gamma;
  auto pbc = PeriodicBoundaries(lengths, angles, false);
  BondOrderCollection bo = Utils::BondDetector::detectBonds(atoms, pbc);
  ASSERT_EQ(bo.getSystemSize(), refBO.getSystemSize());
  for (int i = 0; i < bo.getSystemSize(); ++i) {
    for (int j = 0; j < bo.getSystemSize(); ++j) {
      if (refBO.getOrder(i, j) > 0.1)
        ASSERT_EQ(bo.getOrder(i, j), 1);
      else
        ASSERT_EQ(bo.getOrder(i, j), 0);
    }
  }
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
