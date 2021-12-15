/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/Solvation/SoluteSolventComplex.h"
#include "Utils/Solvation/RandomIndexGenerator.h"
#include <Core/Log.h>

namespace Scine {
namespace Utils {
namespace SoluteSolventComplex {

std::vector<std::vector<AtomCollection>> solvate(const AtomCollection& soluteComplex, int soluteSize,
                                                 const AtomCollection& solvent, int numSolvents, int seed, int resolution,
                                                 double solventOffset, double maxDistance, double stepSize,
                                                 int numRotamers, bool strategicSolv, double coverageThreshold) {
  return SoluteSolventComplex::solvate(soluteComplex, soluteSize, solvent, numSolvents, std::numeric_limits<int>::max(),
                                       seed, resolution, solventOffset, maxDistance, stepSize, numRotamers,
                                       strategicSolv, coverageThreshold);
}

std::vector<std::vector<AtomCollection>> solvateShells(const AtomCollection& soluteComplex, int soluteSize,
                                                       const AtomCollection& solvent, int numShells, int seed, int resolution,
                                                       double solventOffset, double maxDistance, double stepSize,
                                                       int numRotamers, bool strategicSolv, double coverageThreshold) {
  return SoluteSolventComplex::solvate(soluteComplex, soluteSize, solvent, std::numeric_limits<int>::max(), numShells,
                                       seed, resolution, solventOffset, maxDistance, stepSize, numRotamers,
                                       strategicSolv, coverageThreshold);
}

std::vector<std::vector<AtomCollection>>
solvate(const AtomCollection& soluteComplex, const int soluteSize, const AtomCollection& solvent, const int numSolvents,
        const int numShells, const int seed, const int resolution, const double solventOffset, const double maxDistance,
        const double stepSize, const int numRotamers, const bool strategicSolv, const double coverageThreshold) {
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
  // save original number of sites and initialize coverag
  double sizeVisibleComplexSites = visibleComplexSites.size();
  bool resetSize = false;
  double coverage = 0.0;
  // loop till required number of solvent molecules has been added
  while (numAddedSolvents < numSolvents && numAddedShells < numShells) {
    bool addition = false;
    // check, if previous addition round was successful; if yes, recalculate surface
    if (prevAddition && (numAddedSolvents % reCalcTrigger == 0)) {
      visibleComplexSites = MolecularSurface::getVisibleMolecularSurface(complex, soluteStart, soluteEnd, resolution);
      // recalculate number of visible complex sites to calculate coverage
      if (resetSize) {
        sizeVisibleComplexSites = visibleComplexSites.size();
        resetSize = false;
      }
      // calculate coverage of current solute
      coverage = 1.0 - visibleComplexSites.size() / sizeVisibleComplexSites;
      if (strategicSolv) {
        reCalcTrigger = SoluteSolventComplex::solvationStrategy(visibleComplexSites.size());
      }
      // shuffle visible solute sites
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
    if (visibleComplexSites.empty() || coverage >= coverageThreshold) {
      dMin = solventOffset;
      dMax = maxDistance;
      soluteStart = soluteEnd;
      soluteEnd = complex.size();
      resetSize = true;
      prevAddition = true;
      reCalcTrigger = 1;
      coverage = 0.0;
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

AtomCollection mergeAtomCollectionVector(const std::vector<AtomCollection>& atomCollList) {
  AtomCollection merged;
  for (const auto& atoms : atomCollList) {
    merged += atoms;
  }
  return merged;
}

AtomCollection mergeSolventShellVector(const std::vector<std::vector<AtomCollection>>& shellVector) {
  AtomCollection merged;
  for (const auto& atomCollVector : shellVector) {
    merged += SoluteSolventComplex::mergeAtomCollectionVector(atomCollVector);
  }
  return merged;
}

std::vector<int> transferSolventShellVector(const std::vector<std::vector<AtomCollection>>& shellVector) {
  std::vector<int> solventSizeVector;
  for (const auto& shell : shellVector) {
    std::vector<int> tmpSolventSizeVector(shell.size(), 0);
    int count = 0;
    for (const auto& solvent : shell) {
      tmpSolventSizeVector.at(count) = solvent.size();
      count += 1;
    }
    solventSizeVector.insert(solventSizeVector.end(), tmpSolventSizeVector.begin(), tmpSolventSizeVector.end());
  }
  return solventSizeVector;
}

std::vector<std::vector<AtomCollection>>
giveSolventShellVector(const AtomCollection& complex, int soluteSize, const std::vector<int>& solventSizeVector,
                       int resolution, Core::Log& log, bool strategicSolv, double coverageThreshold) {
  int soluteStart = 0;
  int soluteEnd = soluteSize;
  int numAddedSolvents = 0;
  // trigger to force recalculation of the visible surface
  int reCalcTrigger = 1;
  AtomCollection tmpComplex;
  // write solute from complex into tmpSolute
  for (int i = 0; i < soluteSize; i++) {
    tmpComplex.push_back(complex.at(i));
  }
  // set up empty solvent Vector and shell vector
  std::vector<AtomCollection> solventVector;
  std::vector<std::vector<AtomCollection>> shellVector;

  int complexIndex = soluteSize;
  auto visibleComplexSites = MolecularSurface::getVisibleMolecularSurface(tmpComplex, soluteStart, soluteEnd, resolution);
  double sizeVisibleComplexSites = visibleComplexSites.size();
  double coverage = 0.0;
  // loop over the solvent size in solvent size vector
  for (const auto& tmpSize : solventSizeVector) {
    AtomCollection tmpSolvent;
    // get one solvent from complex by current solvent size (tmpSize) and current position (complexIndex)
    while (complexIndex < tmpComplex.size() + tmpSize) {
      tmpSolvent.push_back(complex.at(complexIndex));
      complexIndex += 1;
    }
    // add temporary solvent to current complex
    tmpComplex += tmpSolvent;
    numAddedSolvents += 1;
    // write temporary solvent to solvent vector
    solventVector.push_back(tmpSolvent);
    // recalculate visible surface sites with or without solvation strategy
    if (numAddedSolvents % reCalcTrigger == 0) {
      visibleComplexSites = MolecularSurface::getVisibleMolecularSurface(tmpComplex, soluteStart, soluteEnd, resolution);
      coverage = 1.0 - visibleComplexSites.size() / sizeVisibleComplexSites;
      if (strategicSolv) {
        reCalcTrigger = SoluteSolventComplex::solvationStrategy(visibleComplexSites.size());
      }
    }
    // redefine solute in the case the current solute is solvated
    if ((coverage >= coverageThreshold) || (visibleComplexSites.empty())) {
      soluteStart = soluteEnd;
      soluteEnd = tmpComplex.size();
      visibleComplexSites = MolecularSurface::getVisibleMolecularSurface(tmpComplex, soluteStart, soluteEnd, resolution);
      sizeVisibleComplexSites = visibleComplexSites.size();
      coverage = 0.0;
      shellVector.push_back(solventVector);
      solventVector.clear();
    }
  }
  // add last solvent vector, if it has not been added so far
  if (!visibleComplexSites.empty()) {
    shellVector.push_back(solventVector);
  }
  log.output << "Coverage (Threshold " << coverageThreshold * 100 << " %) of last solute (size = " << soluteEnd
             << "):\t" << coverage << Core::Log::nl;
  return shellVector;
}

int solvationStrategy(int numberSurfPoints) {
  // case if number of surface points is larger than 10000
  if (numberSurfPoints >= 10000) {
    double exponent = log10(numberSurfPoints) - 2;
    return std::pow(10, static_cast<int>(exponent));
  }

  if (numberSurfPoints < 501) {
    // value for 201 - 500
    if (numberSurfPoints > 200) {
      return 25;
    }

    // value for 1 - 20
    if (numberSurfPoints < 21) {
      return 1;
    }

    // value for 21 - 100
    if (numberSurfPoints < 101) {
      return 5;
    }

    // value for 101 - 200
    return 10;
  }

  // value for 501 - 9999
  return 50;
}

bool checkDistances(const AtomCollection& atoms1, const AtomCollection& atoms2) {
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

PositionCollection arrange(const Position& surfPoint1, const Position& surfNormal1, const Position& surfPoint2,
                           const Position& surfNormal2, const PositionCollection& atoms2Position, const double distance) {
  // Final position of second reactive site
  const Position finalPosition = surfPoint1 + surfNormal1 * distance;

  PositionCollection tempPositions(atoms2Position.rows(), 3);
  const Position surfPoint2Rev = -1.0 * surfPoint2;
  // shift surface point of atoms2 to origin
  tempPositions = Geometry::Manipulations::translatePositions(atoms2Position, surfPoint2Rev);
  // rotate postitions to face reactive site of atom 1 at the origin
  const Position pointOfRot(0, 0, 0);
  const Position surfNormal1Rev = -1.0 * surfNormal1;
  tempPositions = Geometry::Manipulations::rotatePositions(tempPositions, surfNormal2, surfNormal1Rev, pointOfRot);
  // shift to final position
  tempPositions = Geometry::Manipulations::translatePositions(tempPositions, finalPosition);

  return tempPositions;
}

bool add(AtomCollection& complex, const AtomCollection& additive, const MolecularSurface::SurfaceSite& complexSurfSite,
         const MolecularSurface::SurfaceSite& additiveSurfSite, double minDistance, double maxDistance,
         double incrementDistance, int numRotationAttempts) {
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
      PositionCollection rotamer =
          Geometry::Manipulations::rotatePositions(additiveArranged, complexSurfNorm, angle, complexSurfPos);
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

} // namespace SoluteSolventComplex
} /* namespace Utils */
} /* namespace Scine */
