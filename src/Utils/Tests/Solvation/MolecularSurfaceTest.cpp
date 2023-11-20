/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Geometry.h"
#include "Utils/IO/ChemicalFileFormats/XyzStreamHandler.h"
#include "Utils/Solvation/MolecularSurface.h"
#include "Utils/Solvation/SurfaceSite.h"
#include "Utils/Typenames.h"
#include <gmock/gmock.h>
#include <fstream>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

TEST(MolecularSurfaceTest, GivesFibonacciInputPoints) {
  int numberOfSurfPoints = 500;
  PositionCollection testSphere = MolecularSurface::FibonacciSphere(numberOfSurfPoints);
  ASSERT_THAT(testSphere.rows(), numberOfSurfPoints);
}

TEST(MolecularSurfaceTest, AreFibonacciPointsOnUnitSphere) {
  int numberOfSurfPoints = 500;
  PositionCollection testSphere = MolecularSurface::FibonacciSphere(numberOfSurfPoints);
  Position center = {0, 0, 0};
  double sumDist = 0;
  for (int i = 0; i < testSphere.rows(); i++) {
    auto dist = (testSphere.row(i) - center).norm();
    sumDist += dist;
  }
  ASSERT_THAT(sumDist, numberOfSurfPoints);
}

TEST(MolecularSurfaceTest, DoesSurfaceSiteConstructDefaultCorrectly) {
  MolecularSurface::SurfaceSite testSurfSite;

  ASSERT_THAT(testSurfSite.position, Position(0, 0, 1));
  ASSERT_THAT(testSurfSite.normal, Position(0, 0, 1));
}

TEST(MolecularSurfaceTest, DoesSurfaceSiteConstructerWork) {
  MolecularSurface::SurfaceSite testSurfaceSite(Position(2, 2, 2), Position(0, 0, 0));

  ASSERT_THAT(testSurfaceSite.position, Position(2, 2, 2));
  ASSERT_THAT(testSurfaceSite.normal.norm(), 1.0);
  ASSERT_THAT(testSurfaceSite.normal, Position(2, 2, 2) / sqrt(12.0));
}

TEST(MolecularSurfaceTest, GivesUnprunedAtomSurfceInputPoints) {
  int numberOfSurfPoints = 500;

  Position atomPosition(0, 0, 0);
  Atom testAtom{ElementType::Xe, atomPosition};

  std::vector<MolecularSurface::SurfaceSite> testAtomSurf =
      MolecularSurface::getUnprunedAtomSurface(testAtom, numberOfSurfPoints);

  ASSERT_THAT(testAtomSurf.size(), numberOfSurfPoints);
}

TEST(MolecularSurfaceTest, IsSurfaceOnVdwRadius) {
  int numberOfSurfPoints = 500;

  Position atomPosition(0, 0, 0);
  Atom testAtom{ElementType::Xe, atomPosition};

  std::vector<MolecularSurface::SurfaceSite> testAtomSurf =
      MolecularSurface::getUnprunedAtomSurface(testAtom, numberOfSurfPoints);

  for (auto& surfSite : testAtomSurf) {
    double distance = (surfSite.position - atomPosition).norm();
    ASSERT_TRUE(distance - ElementInfo::vdwRadius(testAtom.getElementType()) < 1e-12);
  }
}

TEST(MolecularSurfaceTest, PrunesGetAtomBuriedSurfaceSitesCorrect) {
  int surfPoints = 500;

  AtomCollection testMolecule;

  testMolecule.push_back(Atom{ElementType::He, Position(0, 0, 0)});
  testMolecule.push_back(Atom{ElementType::Xe, Position(0, 0, 0)});

  std::vector<MolecularSurface::SurfaceSite> prunedSurfHe;
  std::vector<MolecularSurface::SurfaceSite> prunedSurfXe;

  for (int i = 0; i < testMolecule.size(); i++) {
    auto atom = testMolecule.at(i);
    auto element = atom.getElementType();

    if (element == ElementType::He) {
      prunedSurfHe = MolecularSurface::getPrunedAtomSurface(0, testMolecule, surfPoints);
      ASSERT_THAT(prunedSurfHe.size(), 0);
    }
    if (element == ElementType::Xe) {
      prunedSurfXe = MolecularSurface::getPrunedAtomSurface(1, testMolecule, surfPoints);
      ASSERT_THAT(prunedSurfXe.size(), surfPoints);
    }
  }
}

TEST(MolecularSurfaceTest, GetMolecularSurfaceCorrect) {
  int surfPoints = 64;

  AtomCollection testMolecule;

  for (int i = 0; i < 3; i++) {
    if (i % 2) {
      testMolecule.push_back(Atom{ElementType::Xe, Position(3 * i, 0, 0)});
    }
    else {
      testMolecule.push_back(Atom{ElementType::He, Position(3 * i, 3, 0)});
    }
  }

  std::vector<MolecularSurface::SurfaceSite> testMoleculeSurf =
      MolecularSurface::getPrunedMolecularSurface(testMolecule, surfPoints);

  ASSERT_THAT(testMoleculeSurf.size(), 139);
}

TEST(MolecularSurfaceTest, DoesRayMissSphere) {
  // surface facing sphere, ray should hit sphere
  auto surfSite1 = MolecularSurface::SurfaceSite(Position(2, 2, 2), Position(3, 3, 3));
  auto test1 = MolecularSurface::rayMissesSphere(surfSite1, Position(0, 0, 0), 1);
  // other site of surface, ray should miss sphere
  auto surfSite2 = MolecularSurface::SurfaceSite(Position(2, 2, 2), Position(1.5, 1.5, 1.5));
  auto test2 = MolecularSurface::rayMissesSphere(surfSite2, Position(0, 0, 0), 1);
  // ray should miss sphere, wrong direction
  auto surfSite3 = MolecularSurface::SurfaceSite(Position(2, 2, 4), Position(3, 3, 3));
  auto test3 = MolecularSurface::rayMissesSphere(surfSite3, Position(0, 0, 0), 1);
  // ray should just touch sphere
  auto surfSite4 = MolecularSurface::SurfaceSite(Position(8, 0, 1), Position(9, 0, 1));
  auto test4 = MolecularSurface::rayMissesSphere(surfSite4, Position(0, 0, 0), 1);

  ASSERT_THAT(test1, false);
  ASSERT_THAT(test2, true);
  ASSERT_THAT(test3, true);
  ASSERT_THAT(test4, false);
}

TEST(MolecularSurfaceTest, RayMissSphereIdentifiesOwnSphere) {
  // surface on own sphere, should be true
  auto surfSite1 = MolecularSurface::SurfaceSite(Position(4, 3, 0), Position(0, 0, 0));
  auto test1 = MolecularSurface::rayMissesSphere(surfSite1, Position(0, 0, 0), 5.0);
  // surface touching sphere of same radius as own sphere, should be false
  auto surfSite2 = MolecularSurface::SurfaceSite(Position(4, 3, 0), Position(8, 6, 0));
  auto test2 = MolecularSurface::rayMissesSphere(surfSite2, Position(0, 0, 0), 5.0);

  ASSERT_THAT(test1, true);
  ASSERT_THAT(test2, false);
}

TEST(MolecularSurfaceTest, GetVisibleMolecularSurfaceCorrect) {
  int surfPoints = 32;
  int k = 0;
  AtomCollection testMolecule;
  for (int i = 0; i < 3; i++) {
    if (i % 2) {
      testMolecule.push_back(Atom{ElementType::Xe, Position(0 * i, 0, 0)});
    }
    else {
      testMolecule.push_back(Atom{ElementType::F, Position(1.5 * 2 * std::pow(-1, k), 1.5 * 2, 0)});
      k += 1;
    }
  }

  PositionCollection solventPoints = MolecularSurface::FibonacciSphere(30);

  for (int i = 0; i < solventPoints.rows(); i++) {
    auto solventPos = solventPoints.row(i) * 10.0;
    testMolecule.push_back(Atom{ElementType::He, solventPos});
  }

  auto visibleSitesTestMolecule = MolecularSurface::getVisibleMolecularSurface(testMolecule, 0, 3, surfPoints);

  ASSERT_THAT(visibleSitesTestMolecule.size(), 27);
}

TEST(MolecularSurfaceTest, GetVisibleMolecularSurfaceSoluteSeeingEachOther) {
  int surfPoints = 256;
  // Test case set up that solute atoms 'see' each other, meaning they are within the distance threshold of 10 bohr
  double xCoord = 7.6;

  AtomCollection testMolecule;
  testMolecule.push_back(Atom{ElementType::He, Position{-xCoord, 0, 0}});
  testMolecule.push_back(Atom{ElementType::He, Position{xCoord, 0, 0}});

  auto visibleInternalSoluteSites =
      MolecularSurface::getVisibleMolecularSurface(testMolecule, 0, testMolecule.size(), surfPoints);

  testMolecule.push_back(Atom{ElementType::Cs, Position{0, 0, 0}});

  auto visibleSoluteSites = MolecularSurface::getVisibleMolecularSurface(testMolecule, 0, 2, surfPoints);

  ASSERT_THAT(static_cast<int>(visibleInternalSoluteSites.size()) < surfPoints * 2, true);
  ASSERT_THAT(visibleSoluteSites.size() < visibleInternalSoluteSites.size(), true);
}

TEST(MolecularSurfaceTest, GetVisibleMolecularSurfaceSoluteNotSeeingEachOther) {
  int surfPoints = 256;
  // Test case set up that solute atoms do not 'see' each other, meaning they are just outside the distance threshold of
  // 10 bohr
  double xCoord = 7.65;

  AtomCollection testMolecule;
  testMolecule.push_back(Atom{ElementType::He, Position{-xCoord, 0, 0}});
  testMolecule.push_back(Atom{ElementType::He, Position{xCoord, 0, 0}});

  auto visibleInternalSoluteSites =
      MolecularSurface::getVisibleMolecularSurface(testMolecule, 0, testMolecule.size(), surfPoints);

  testMolecule.push_back(Atom{ElementType::Cs, Position{0, 0, 0}});

  auto visibleSoluteSites = MolecularSurface::getVisibleMolecularSurface(testMolecule, 0, 2, surfPoints);

  ASSERT_THAT(visibleInternalSoluteSites.size(), 2 * surfPoints);
  ASSERT_THAT(visibleSoluteSites.size() < visibleInternalSoluteSites.size(), true);
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */

// Code to Export molecule and its surface to files
// std::ofstream fout2;
// fout2.open("/home/cds/tpaul/testMolecule.xyz");
// XyzStreamHandler::write(fout2, testMolecule);
// fout2.close();
//
// std::ofstream fout;
// fout.open("/home/cds/tpaul/testSurface.xyz");
// MolecularSurface::writeSurface(fout, testAtomSurfPos);
// fout.close();
//
// std::ofstream fout1;
// fout1.open("/home/cds/tpaul/testSurfaceNorm.xyz");
// MolecularSurface::writeSurface(fout1, testAtomSurfNorm);
// fout1.close();
