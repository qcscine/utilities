/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/Solvation/SoluteSolventComplex.h"
#include "Utils/Solvation/RandomIndexGenerator.h"

namespace Scine {
namespace Utils {

std::vector<std::vector<AtomCollection>>
SoluteSolventComplex::solvate(const AtomCollection& soluteComplex, int soluteSize, const AtomCollection& solvent,
                              int numSolvents, int seed, int resolution, double solventOffset, double maxDistance,
                              double stepSize, int numRotamers, bool strategicSolv) {
  return SoluteSolventComplex::solvate(soluteComplex, soluteSize, solvent, numSolvents, std::numeric_limits<int>::max(),
                                       seed, resolution, solventOffset, maxDistance, stepSize, numRotamers, strategicSolv);
}

std::vector<std::vector<AtomCollection>>
SoluteSolventComplex::solvateShells(const AtomCollection& soluteComplex, int soluteSize, const AtomCollection& solvent,
                                    int numShells, int seed, int resolution, double solventOffset, double maxDistance,
                                    double stepSize, int numRotamers, bool strategicSolv) {
  return SoluteSolventComplex::solvate(soluteComplex, soluteSize, solvent, std::numeric_limits<int>::max(), numShells,
                                       seed, resolution, solventOffset, maxDistance, stepSize, numRotamers, strategicSolv);
}

std::vector<std::vector<AtomCollection>>
SoluteSolventComplex::solvate(const AtomCollection& soluteComplex, int soluteSize, const AtomCollection& solvent,
                              int numSolvents, int numShells, int seed, int resolution, double solventOffset,
                              double maxDistance, double stepSize, int numRotamers, bool strategicSolv) {
  auto solventSites = MolecularSurface::getPrunedMolecularSurface(solvent, resolution);

  double dMin = solventOffset;
  double dMax = maxDistance;

  auto complex = soluteComplex;
  std::vector<AtomCollection> solventVector;
  std::vector<std::vector<AtomCollection>> shellVector;
  int numAddedSolvents = 0;
  int numAddedShells = 0;
  int soluteStart = 0;
  int soluteEnd = soluteSize;
  // random number generator with given seed
  std::mt19937 shuffleGen(seed);
  // set up index generator with seed and solvent
  SoluteSolventComplex::RandomIndexGenerator randSolventSite(solventSites.size(), seed);

  // bool to remember if previous addition was successful or not
  bool prevAddition = false;
  int reCalcCounter = 0;
  int failedAdd = 0;
  // trigger to force recalculation of the visible surface
  int reCalcTrigger = 1;
  // calculate first surface before adding anything
  auto visibleComplexSites = MolecularSurface::getVisibleMolecularSurface(complex, soluteStart, soluteEnd, resolution);

  // loop till required number of solvent molecules has been added
  while (numAddedSolvents < numSolvents && numAddedShells < numShells) {
    bool addition = false;
    // check, if previous addition round was successful; if yes, recalculate surface
    if (prevAddition && (numAddedSolvents % reCalcTrigger == 0)) {
      visibleComplexSites = MolecularSurface::getVisibleMolecularSurface(complex, soluteStart, soluteEnd, resolution);
      if (strategicSolv) {
        reCalcTrigger = SoluteSolventComplex::solvationStrategy(visibleComplexSites.size());
      }
      std::shuffle(std::begin(visibleComplexSites), std::end(visibleComplexSites), shuffleGen);
      reCalcCounter += 1;
    }
    // choose random solvent surface index
    int randSolventSiteIndex = randSolventSite.next();
    // loop over visible complex sites and try to add solvent molecule
    for (const auto& visibleComplexSite : visibleComplexSites) {
      addition = SoluteSolventComplex::add(complex, solvent, visibleComplexSite, solventSites.at(randSolventSiteIndex),
                                           dMin, dMax, stepSize, numRotamers);
      // if addition was successful, break loop over complex sites and try to add next solvent molecule
      if (addition) {
        numAddedSolvents += 1;
        prevAddition = true;
        AtomCollection tmpSolvent;
        int tmpIndex = 0;
        // extract lastly added solvent and push to solvent list
        for (int i = complex.size() - solvent.size(); i < complex.size(); i++) {
          tmpSolvent.push_back(complex.at(i));
          tmpIndex += 1;
        }
        solventVector.push_back(tmpSolvent);
        break;
      }
    }
    // if addition over all sites was not successful, try next round by increasing the dMin and dMax values
    if (!addition) {
      dMin = dMax;
      dMax += maxDistance;
      prevAddition = false;
      failedAdd += 1;
    }
    // if no more complex sites are visible, return solventVector
    if (visibleComplexSites.empty()) {
      dMin = solventOffset;
      dMax = maxDistance;
      soluteStart = soluteEnd;
      soluteEnd = complex.size();
      prevAddition = true;
      numAddedShells += 1;
      shellVector.push_back(solventVector);
      solventVector.clear();
    }
  }
  // only add last solventVector, if it has not been added due to a full shell
  if (!visibleComplexSites.empty()) {
    shellVector.push_back(solventVector);
  }
  return shellVector;
}

AtomCollection SoluteSolventComplex::mergeAtomCollectionVector(const std::vector<AtomCollection>& atomCollList) {
  AtomCollection merged;
  for (const auto& atoms : atomCollList) {
    merged += atoms;
  }
  return merged;
}

AtomCollection SoluteSolventComplex::mergeSolventShellVector(const std::vector<std::vector<AtomCollection>>& shellVector) {
  AtomCollection merged;
  for (const auto& atomCollVector : shellVector) {
    merged += SoluteSolventComplex::mergeAtomCollectionVector(atomCollVector);
  }
  return merged;
}

int SoluteSolventComplex::solvationStrategy(int numberSurfPoints) {
  // case if number of surface points is larger than 10000
  if (numberSurfPoints >= 10000) {
    double exponent = log10(numberSurfPoints) - 2;
    return pow(10, (int)exponent);
  }
  else {
    if (numberSurfPoints < 501) {
      // value for 201 - 500
      if (numberSurfPoints > 200) {
        return 25;
      }
      else {
        // value for 1 - 20
        if (numberSurfPoints < 21) {
          return 1;
        }
        // value for 21 - 100
        if (numberSurfPoints < 101) {
          return 5;
        }
        // value for 101 - 200
        else {
          return 10;
        }
      }
    }
    // value for 501 - 9999
    else {
      return 50;
    }
  }
}

bool SoluteSolventComplex::checkDistances(const AtomCollection& atoms1, const AtomCollection& atoms2) {
  // loop over atom collection, assuming atoms2 is smaller than atoms1
  for (const auto& atom2 : atoms2) {
    const Position& posAtom2 = atom2.getPosition();
    double r2 = ElementInfo::vdwRadius(atom2.getElementType());
    // loop over atom collection atoms1
    for (const auto& atom1 : atoms1) {
      double distance = (atom1.getPosition() - posAtom2).norm();
      // check, if the atoms are close; threshold set to two times the approx. radius of Cs
      if ((distance - 13.0) < 1e-12) {
        double r1 = ElementInfo::vdwRadius(atom1.getElementType());
        // if they are close, check if they overlap
        if (distance < r1 + r2) {
          return false;
        }
      }
    }
  }
  return true;
}

PositionCollection SoluteSolventComplex::arrange(Position surfPoint1, Position surfNormal1, Position surfPoint2,
                                                 Position surfNormal2, const PositionCollection& atoms2Position,
                                                 double distance) {
  // Final position of second reactive site
  Position finalPosition = surfPoint1 + surfNormal1 * distance;

  PositionCollection tempPositions(atoms2Position.rows(), 3);
  Position surfPoint2Rev = -1.0 * surfPoint2;
  // shift surface point of atoms2 to origin
  tempPositions = Geometry::translatePositions(atoms2Position, surfPoint2Rev);
  // rotate postitions to face reactive site of atom 1 at the origin
  Position pointOfRot(0, 0, 0);
  Position surfNormal1Rev = -1.0 * surfNormal1;
  tempPositions = Geometry::rotatePositions(tempPositions, surfNormal2, surfNormal1Rev, pointOfRot);
  // shift to final position
  tempPositions = Geometry::translatePositions(tempPositions, finalPosition);

  return tempPositions;
}

bool SoluteSolventComplex::add(AtomCollection& complex, const AtomCollection& additive,
                               const MolecularSurface::SurfaceSite& complexSurfSite,
                               const MolecularSurface::SurfaceSite& additiveSurfSite, double minDistance,
                               double maxDistance, double incrementDistance, int numRotationAttempts) {
  // local copy of surface position and surface normal vector
  auto complexSurfPos = complexSurfSite.position;
  auto additiveSurfPos = additiveSurfSite.position;
  auto complexSurfNorm = complexSurfSite.normal;
  auto additiveSurfNorm = additiveSurfSite.normal;

  double distance = minDistance;
  // loop until as long as distance is equal or smaller than maxDistance
  while (distance <= maxDistance) {
    // trial arrangement of additive relative to surface site of complex
    PositionCollection additiveArranged =
        arrange(complexSurfPos, complexSurfNorm, additiveSurfPos, additiveSurfNorm, additive.getPositions(), distance);
    // loop over number of rotations
    for (int rotation = 0; rotation < numRotationAttempts; rotation++) {
      double angle = 2 * Constants::pi * rotation / numRotationAttempts;
      // rotated trial arrangement of additive
      PositionCollection rotamer = Geometry::rotatePositions(additiveArranged, complexSurfNorm, angle, complexSurfPos);
      // set AtomCollection of additive together
      AtomCollection trialAtoms(additive.getElements(), rotamer);
      // check, if trial atoms overlap with any atoms of complex; if true, none overlap and the two atom collections are
      // combined
      if (checkDistances(complex, trialAtoms)) {
        complex += trialAtoms;
        return true;
      }
    }

    distance += incrementDistance;
  }

  return false;
}

} /* namespace Utils */
} /* namespace Scine */
