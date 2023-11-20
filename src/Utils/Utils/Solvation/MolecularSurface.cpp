/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Solvation/MolecularSurface.h"
#include "Utils/IO/ChemicalFileFormats/XyzStreamHandler.h"

namespace Scine {
namespace Utils {

PositionCollection MolecularSurface::FibonacciSphere(int numberOfPoints) {
  auto offset = 2.0 / numberOfPoints;
  auto increment = Constants::pi * (3.0 - std::sqrt(5.0));

  PositionCollection fiboSurfacePoints(numberOfPoints, 3);

  for (int i = 0; i < numberOfPoints; i++) {
    auto u = ((i * offset) - 1) + (offset / 2);
    auto r = std::sqrt(1 - u * u);
    auto phi = ((i + 1) % numberOfPoints) * increment;

    fiboSurfacePoints.row(i) = Position(std::cos(phi) * r, u * 1, std::sin(phi) * r);
  }
  return fiboSurfacePoints;
}

std::vector<MolecularSurface::SurfaceSite> MolecularSurface::getUnprunedAtomSurface(const Atom& atom, int atomSurfPoints) {
  // set up fibonacci sphere of given resolution
  auto initialPoints = FibonacciSphere(atomSurfPoints);
  // copy information of atom
  const auto& atomPosition = atom.getPosition();
  auto element = atom.getElementType();
  auto r = ElementInfo::vdwRadius(element);

  std::vector<MolecularSurface::SurfaceSite> unprunedAtomSurface(atomSurfPoints);

  // build unpruned surface of atom
  for (int j = 0; j < initialPoints.rows(); j++) {
    MolecularSurface::SurfaceSite siteJ(initialPoints.row(j) * r + atomPosition, atomPosition);
    unprunedAtomSurface.at(j) = siteJ;
  }

  return unprunedAtomSurface;
}

std::vector<MolecularSurface::SurfaceSite> MolecularSurface::getPrunedAtomSurface(int atomIndex, const AtomCollection& atoms,
                                                                                  int atomSurfPoints) {
  //  write indices of close atoms to list
  std::vector<int> closeAtoms;
  // position of atom given by index
  Position atomIndexPosition = atoms.at(atomIndex).getPosition();
  // identify atoms which are within a sphere of 10
  for (int i = 0; i < atoms.size(); i++) {
    double distance = (atoms.at(i).getPosition() - atomIndexPosition).norm();
    if ((distance - 10.0) < 1e-12 && atomIndex != i) {
      closeAtoms.push_back(i);
    }
  }
  //  calculate unpruned surface points of atom with index atomIndex
  std::vector<MolecularSurface::SurfaceSite> unprunedAtomPoints = getUnprunedAtomSurface(atoms.at(atomIndex), atomSurfPoints);
  // bool vector to account for remaining atom points
  std::vector<bool> remainingAtomPointsBool(unprunedAtomPoints.size(), true);
  int remainingAtomPointsSize = 0;
  // loop over all atom points
  for (int j = 0; j < static_cast<int>(unprunedAtomPoints.size()); j++) {
    bool noHit = true;
    // loop over the atoms close to given atom
    for (const auto& closeAtom : closeAtoms) {
      double distance = (atoms.at(closeAtom).getPosition() - unprunedAtomPoints.at(j).position).norm();
      auto element = atoms.at(closeAtom).getElementType();
      auto vdw = ElementInfo::vdwRadius(element);
      if ((distance - vdw) < 1e-12) {
        noHit = false;
        break;
      }
    }
    // save bool to vector
    remainingAtomPointsBool.at(j) = noHit;
    // increase remaining point counter by 1
    if (noHit) {
      remainingAtomPointsSize += 1;
    }
  }
  // setup vector for remaining atom points
  std::vector<MolecularSurface::SurfaceSite> remainingAtomPoints(remainingAtomPointsSize);
  int indexRemainingAtomsPoint = 0;
  // loop over bool vector to write remaining positions from unpruned atom points to new vector
  for (int i = 0; i < static_cast<int>(remainingAtomPointsBool.size()); i++) {
    if (remainingAtomPointsBool.at(i)) {
      remainingAtomPoints.at(indexRemainingAtomsPoint) = unprunedAtomPoints.at(i);
      indexRemainingAtomsPoint += 1;
    }
  }

  return remainingAtomPoints;
}

std::vector<MolecularSurface::SurfaceSite> MolecularSurface::getPrunedMolecularSurface(const AtomCollection& atoms,
                                                                                       int atomSurfPoints) {
  std::vector<MolecularSurface::SurfaceSite> molecularSurface;
  // loop over atoms in atom collection
  for (int atomIndex = 0; atomIndex < atoms.size(); atomIndex++) {
    // prune atom surface
    std::vector<MolecularSurface::SurfaceSite> prunedAtomSurface = getPrunedAtomSurface(atomIndex, atoms, atomSurfPoints);
    // add to molecular surface
    for (const auto& atomSurfaceSite : prunedAtomSurface) {
      molecularSurface.push_back(atomSurfaceSite);
    }
  }

  return molecularSurface;
}

bool MolecularSurface::rayMissesSphere(const MolecularSurface::SurfaceSite& surfaceSite, const Position& sphereOrigin,
                                       double sphereRad) {
  // calculate vector from sphere origin to surface position
  Position shifted = surfaceSite.position - sphereOrigin;
  // check if surface position and sphere belong to the same atom; if so, ray misses sphere is true
  if ((shifted / sphereRad - surfaceSite.normal).norm() < 1e-12) {
    return true;
  }
  double b = 2 * surfaceSite.normal.dot(shifted);
  double c = shifted.squaredNorm() - sphereRad * sphereRad;

  double discriminant = b * b - 4 * c;
  if (discriminant < 0.0) {
    return true;
  }
  // only positive solution is acceptable, therefore discriminant has to be smaller than b
  return std::sqrt(discriminant) < b;
}

std::vector<MolecularSurface::SurfaceSite>
MolecularSurface::getVisibleMolecularSurface(const AtomCollection& molecule, int startIndex, int endIndex, int resolution) {
  AtomCollection solute;
  // loop over complex to extract solute as local copy
  for (int i = startIndex; i < endIndex; i++) {
    solute.push_back(molecule.at(i));
  }
  std::vector<MolecularSurface::SurfaceSite> visibleMolecularSurface;
  // loop over atoms of solute
  for (int j = 0; j < solute.size(); j++) {
    std::vector<MolecularSurface::SurfaceSite> prunedAtomSites = getPrunedAtomSurface(j, solute, resolution);

    // bool for visible sites of solute atom j
    std::vector<bool> visibleAtomSitesBool(prunedAtomSites.size(), true);
    // loop over pruned atom sites
    for (int k = 0; k < static_cast<int>(prunedAtomSites.size()); k++) {
      auto prunedAtomSite = prunedAtomSites.at(k);
      // loop over all atoms in complex to check if pruned site hits them
      for (int l = 0; l < molecule.size(); l++) {
        // check if l is worth calculating rayMisses function
        double vdwSphere = ElementInfo::vdwRadius(molecule.at(l).getElementType());
        // only perform ray calculation if atom l is part of solvent OR the distance between the surface of solute atom
        // l and the surface point k is smaller than 10 bohr
        if ((l >= endIndex) ||
            (l < endIndex && (molecule.at(l).getPosition() - prunedAtomSite.position).norm() - vdwSphere < 10.0)) {
          // check, if ray from surface point k hits sphere of atom
          bool rayMisses = rayMissesSphere(prunedAtomSite, molecule.at(l).getPosition(), vdwSphere);
          // break loop if ray hits a sphere;
          if (!rayMisses) {
            visibleAtomSitesBool.at(k) = rayMisses;
            break;
          }
        }
      }
    }
    // write visible atom sites to position collection of molecule
    for (int m = 0; m < static_cast<int>(visibleAtomSitesBool.size()); m++) {
      if (visibleAtomSitesBool.at(m)) {
        visibleMolecularSurface.push_back(prunedAtomSites.at(m));
      }
    }
  }

  return visibleMolecularSurface;
}

void MolecularSurface::writeSurface(std::ostream& os, std::vector<MolecularSurface::SurfaceSite> surface) {
  AtomCollection surfaceInHAtoms(surface.size());

  for (int i = 0; i < static_cast<int>(surface.size()); i++) {
    surfaceInHAtoms.setElement(i, ElementType::H);
    surfaceInHAtoms.setPosition(i, surface.at(i).position);
  }
  XyzStreamHandler::write(os, surfaceInHAtoms);
}

} /* namespace Utils */
} /* namespace Scine */
