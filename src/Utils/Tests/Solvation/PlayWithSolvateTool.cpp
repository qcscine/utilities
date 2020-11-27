/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
//#include "Utils/IO/ChemicalFileFormats/XyzStreamHandler.h"
//#include "Utils/Solvation/SoluteSolventComplex.h"
//#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
//#include <gmock/gmock.h>
//#include <fstream>
//#include <iostream>
//#include <random>
//
// using namespace testing;
// namespace Scine {
// namespace Utils {
// namespace Tests {
//
//// Better compile in release mode
//
// TEST(SoluteSolventComplexTest, DoesSolventListWork) {
//  ElementTypeCollection waterEC(3);
//  waterEC.at(0) = ElementType::O;
//  waterEC.at(1) = ElementType::H;
//  waterEC.at(2) = ElementType::H;
//
//  PositionCollection waterPC0(3, 3);
//  waterPC0.row(0) = Position(0, 0, 0);
//  waterPC0.row(1) = Position(-1.63, 1.25, 0);
//  waterPC0.row(2) = Position(1.63, 1.25, 0);
//
//  AtomCollection water(waterEC, waterPC0);
//
//  AtomCollection solute;
//  solute.push_back(Atom(ElementType::Na, Position(0, 0, 0)));
//
//  int sampleRes = 12;
//  int numSolvents = 10000;

//  auto solventList =
//    SoluteSolventComplex::solvateList(solute, solute.size(), water, numSolvents, 4, 5, sampleRes, 0, 7, 1, 3, true);

//  std::ofstream fout2;
//
//  for(auto & i : solventList) {
//    std::string path = "/home/cds/tpaul/structures_test/test_solventList_";
//    path += std::to_string(i.size()) + ".xyz";
//    fout2.open(path);
//    std::cout << i.size() << std::endl;
//    XyzStreamHandler::write(fout2, SoluteSolventComplex::mergeAtomCollectionVector(i));
//    fout2.close();
//  }

//  auto complex = SoluteSolventComplex::solvateShells(solute, solute.size(), water, 4, 5, sampleRes, 0, 7, 1, 3, true);

//  std::string path = "/home/cds/tpaul/structures_test/test_classicComplex";
//  path += ".xyz";
//  fout2.open(path);
//  XyzStreamHandler::write(fout2, complex);
//  fout2.close();

//  std::cout << solventList.size() << std::endl;
//}
//
// TEST(PlayWithSolvateToolTest, PlayWithMethanal) {
//  PositionCollection methanalPC(4, 3);
//  methanalPC.row(0) = Position(2, 0, 0);
//  methanalPC.row(1) = Position(0, 0, 0);
//  methanalPC.row(2) = Position(-1.8, 1.5, 0);
//  methanalPC.row(3) = Position(-1.8, -1.5, 0);
//
//  ElementTypeCollection methanalEC(4);
//  methanalEC.at(0) = ElementType::O;
//  methanalEC.at(1) = ElementType::C;
//  methanalEC.at(2) = ElementType::H;
//  methanalEC.at(3) = ElementType::H;
//
//  AtomCollection methanal(methanalEC, methanalPC);
//
//  const AtomCollection& orgMethanal = methanal;
//
//  ElementTypeCollection waterEC(3);
//  waterEC.at(0) = ElementType::O;
//  waterEC.at(1) = ElementType::H;
//  waterEC.at(2) = ElementType::H;
//
//  PositionCollection waterPC0(3, 3);
//  waterPC0.row(0) = Position(0, 0, 0);
//  waterPC0.row(1) = Position(-1.8, 1.25, 0);
//  waterPC0.row(2) = Position(-1.8, -1.25, 0);
//
//  AtomCollection water(waterEC, waterPC0);
//
//  AtomCollection solventHe;
//  solventHe.push_back(Atom(ElementType::He, Position(0, 0, 0)));
//
//  auto solventComplex_start = SoluteSolventComplex::solvate(methanal, methanal.size(), water, 1, 42, 24, 0, 7, 3, 3,
//  true); auto soluteComplex_start = methanal +
//  SoluteSolventComplex::mergeAtomCollectionVector(solventComplex_start.at(0));
//
//  std::ofstream fout2;
//  for (int i = 0; i < 100; i++) {
//    double minDist = i * 0.1;
//    auto tempSoluteComplex = SoluteSolventComplex::solvate(soluteComplex_start, methanal.size(), water, 1, 42, 24,
//                                                           minDist, minDist + 7, 3, 1, true);
//    std::string path = "/home/cds/tpaul/20FS_ETHZ/SolvationProject/cluster_SmallTestSet_Distance/methanalIn2Water_";
//    path += std::to_string(i) + ".xyz";
//    fout2.open(path);
//    XyzStreamHandler::write(fout2, SoluteSolventComplex::mergeAtomCollectionVector(tempSoluteComplex.at(0)));
//    fout2.close();
//  }
//  auto soluteComplex_18 = SoluteSolventComplex::solvate(methanal, methanal.size(), water, 18, 42, 24, 0, 7, 3, 3,
//  true); for (int i = 24; i < 48; i++) {
//    auto tempSoluteComplex =
//        SoluteSolventComplex::solvate(soluteComplex_18, methanal.size(), water, 2, 42 + i - 24, 24, 0, 7, 3, 3,
//        true);
//    std::string path = "/home/cds/tpaul/20FS_ETHZ/SolvationProject/clusterTestSet/methanalIn20Water_";
//    path += std::to_string(i) + ".xyz";
//    fout2.open(path);
//    XyzStreamHandler::write(fout2, tempSoluteComplex);
//    fout2.close();
//  }
//  auto soluteComplex_17 = SoluteSolventComplex::solvate(methanal, methanal.size(), water, 17, 42, 24, 0, 7, 3, 3,
//  true); for (int i = 48; i < 72; i++) {
//    auto tempSoluteComplex =
//        SoluteSolventComplex::solvate(soluteComplex_17, methanal.size(), water, 3, 42 + i - 48, 24, 0, 7, 3, 3,
//        true);
//    std::string path = "/home/cds/tpaul/20FS_ETHZ/SolvationProject/clusterTestSet/methanalIn20Water_";
//    path += std::to_string(i) + ".xyz";
//    fout2.open(path);
//    XyzStreamHandler::write(fout2, tempSoluteComplex);
//    fout2.close();
//  }
//  auto soluteComplex_16 = SoluteSolventComplex::solvate(methanal, methanal.size(), water, 16, 42, 24, 0, 7, 3, 3,
//  true); for (int i = 72; i < 96; i++) {
//    auto tempSoluteComplex =
//      SoluteSolventComplex::solvate(soluteComplex_16, methanal.size(), water, 4, 42 + i - 72, 24, 0, 7, 3, 3, true);
//    std::string path = "/home/cds/tpaul/20FS_ETHZ/SolvationProject/clusterTestSet/methanalIn20Water_";
//    path += std::to_string(i) + ".xyz";
//    fout2.open(path);
//    XyzStreamHandler::write(fout2, tempSoluteComplex);
//    fout2.close();
//  }
//  auto soluteComplex_15 = SoluteSolventComplex::solvate(methanal, methanal.size(), water, 15, 42, 24, 0, 7, 3, 3,
//  true); for (int i = 96; i < 120; i++) {
//    auto tempSoluteComplex =
//      SoluteSolventComplex::solvate(soluteComplex_15, methanal.size(), water, 5, 42 + i - 96, 24, 0, 7, 3, 3, true);
//    std::string path = "/home/cds/tpaul/20FS_ETHZ/SolvationProject/clusterTestSet/methanalIn20Water_";
//    path += std::to_string(i) + ".xyz";
//    fout2.open(path);
//    XyzStreamHandler::write(fout2, tempSoluteComplex);
//    fout2.close();
//  }
//  for (int i = 120; i < 144; i++) {
//    auto tempSoluteComplex =
//      SoluteSolventComplex::solvate(methanal, methanal.size(), water, 20, 42 + i - 120, 24, 0, 7, 3, 3, true);
//    std::string path = "/home/cds/tpaul/20FS_ETHZ/SolvationProject/clusterTestSet/methanalIn20Water_";
//    path += std::to_string(i) + ".xyz";
//    fout2.open(path);
//    XyzStreamHandler::write(fout2, tempSoluteComplex);
//    fout2.close();
//  }

//  auto soluteComplex2 = SoluteSolventComplex::solvate(methanal, methanal.size(), methanal, 100, 5, 32, 0, 7, 3, 3);

//  auto stillVisibleSites = MolecularSurface::getVisibleMolecularSurface(soluteComplex2, 0, orgMethanal.size(), 32);
//  auto clusterVisibleSites =
//      MolecularSurface::getVisibleMolecularSurface(soluteComplex, orgMethanal.size(), soluteComplex.size(), 32);
//  auto clusterVisibleSites2 = MolecularSurface::getVisibleMolecularSurface(soluteComplex2, orgMethanal.size(), 34 *
//  4, 32);
//
//  std::cout << "Org. Solute:\t" << stillVisibleSites.size() << std::endl;
//  std::cout << "New Complex:\t" << clusterVisibleSites.size() << std::endl;
//  std::cout << "First Shell Complex:\t" << clusterVisibleSites2.size() << std::endl;
//
//
//  fout2.open("/home/cds/tpaul/testSolvate.xyz");
//  XyzStreamHandler::write(fout2, soluteComplex);
//  XyzStreamHandler::write(fout2, soluteComplex2);
//  fout2.close();
//
//  std::ofstream fout;
//  fout.open("/home/cds/tpaul/testSurfaceStillVisible.xyz");
//  MolecularSurface::writeSurface(fout, stillVisibleSites);
//  MolecularSurface::writeSurface(fout, clusterVisibleSites2);
//  fout.close();
//}
//
// TEST(SoluteSolventComplexTest, PlayWithProtein) {
//  auto proteinData = ChemicalFileHandler::read("/home/cds/tpaul/projects/solvation/1L2Y.xyz");
//  auto proteinAtoms = proteinData.first;
//
//  auto benzeneData = ChemicalFileHandler::read("/home/cds/tpaul/benzene.xyz");
//  auto benzeneAtoms = benzeneData.first;
//
//  int sampleRes = 10;
//
//  auto proteinVisibleSurfaceSite =
//      MolecularSurface::getVisibleMolecularSurface(proteinAtoms, 0, proteinAtoms.size(), sampleRes);
//  std::cout << "\tNumber of sites: " << proteinVisibleSurfaceSite.size()
//            << "\t Sites/res: " << proteinVisibleSurfaceSite.size() << std::endl;
//
//  ElementTypeCollection waterEC(3);
//  waterEC.at(0) = ElementType::O;
//  waterEC.at(1) = ElementType::H;
//  waterEC.at(2) = ElementType::H;
//
//  PositionCollection waterPC0(3, 3);
//  waterPC0.row(0) = Position(0, 0, 0);
//  waterPC0.row(1) = Position(-1.63, 1.25, 0);
//  waterPC0.row(2) = Position(1.63, 1.25, 0);
//
//  AtomCollection water(waterEC, waterPC0);
//
//  auto proteinSolventComplex200 = SoluteSolventComplex::solvateShells(proteinAtoms, proteinAtoms.size(), benzeneAtoms,
//                                                                      3, 5, sampleRes, 0, 7, 3, 6, true);
//
//  proteinVisibleSurfaceSite = MolecularSurface::getVisibleMolecularSurface(proteinSolventComplex200,
//  proteinAtoms.size(),
//                                                                           proteinSolventComplex200.size(), sampleRes);
//  std::cout << "Added Molecules: " << (proteinSolventComplex200.size() - proteinAtoms.size()) / 12 << std::endl;
//
//  std::ofstream fout;
//  fout.open("/home/cds/tpaul/testProtein1stShellSurface.xyz");
//  MolecularSurface::writeSurface(fout, proteinVisibleSurfaceSite);
//  fout.close();
//
//  std::ofstream fout2;
//  fout2.open("/home/cds/tpaul/testProteinSolvated.xyz");
//  XyzStreamHandler::write(fout2, proteinSolventComplex200);
//  fout2.close();
//}
//
//} /* namespace Tests */
//} /* namespace Utils */
//} /* namespace Scine */
