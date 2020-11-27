/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILSOS_MOLECULARSURFACE_H
#define UTILSOS_MOLECULARSURFACE_H

#include "Utils/Geometry.h"
#include "Utils/Solvation/SurfaceSite.h"
#include "Utils/Typenames.h"

namespace Scine {
namespace Utils {
/**
 * @brief A approximated molecular surface
 *
 * Namespace to calculate the approximated surface points of a molecule using Fibonacci spheres
 */
namespace MolecularSurface {
/**
 * @brief Build a Fibonacci sphere of N points.
 * @param numberOfPoints The number of points N for the Fibonacci sphere.
 * @return Position Collection with N rows and the coordinates of the N points to build up the sphere.
 */
PositionCollection FibonacciSphere(int numberOfPoints);
/**
 * @brief Build unpruned atom surface around an atom
 * @param atom Atom of interest.
 * @param atomSurfPoints Number of surface points.
 * @return Vector of unpruned surface sites.
 */
std::vector<MolecularSurface::SurfaceSite> getUnprunedAtomSurface(const Atom& atom, int atomSurfPoints);
/**
 * @brief Prune the surface points of one atom.
 * @param atomIndex Index of atom in given AtomCollection.
 * @param atoms Given AtomCollection, corresponding to a molecule.
 * @param atomSurfPoints Surface points per atom, corresponds to resolution.
 * @return Vector of the pruned surface sites of the atom.
 */
std::vector<MolecularSurface::SurfaceSite> getPrunedAtomSurface(int atomIndex, const AtomCollection& atoms, int atomSurfPoints);
/**
 * @brief Gives pruned approximated molecular surface.
 * @param atoms AtomCollection of given molecule of interest.
 * @param atomSurfPoints Unpruned number of surface points per atom, resolution of surface.
 * @return Vector of surface sites of pruned molecular surface.
 */
std::vector<MolecularSurface::SurfaceSite> getPrunedMolecularSurface(const AtomCollection& atoms, int atomSurfPoints);
/**
 * @brief Checks if surface point misses given sphere
 *
 * Based on quadratic equation of (surface.position + t * surface.position.norm - sphere.origin)^2 - sphere.radius^2 = 0,
 * where one solves for t. If t is a real, positive(!) number (including 0), the ray hits a sphere.
 * This can be solved by checking the discriminant.
 *
 * @param surfaceSite Surface site of interest.
 * @param sphereOrigin Position of origin of sphere.
 * @param sphereRad Radius of sphere.
 * @return Bool, true if ray misses, false if ray hits sphere;
 */
bool rayMissesSphere(const MolecularSurface::SurfaceSite& surfaceSite, const Position& sphereOrigin, double sphereRad);
/**
 * @brief Finds surface sites of molecule in complex which do not 'see' other atoms.
 *
 * @param molecule Atom collection of several molecules.
 * @param startIndex Start index of solute in complex.
 * @param endIndex End index of solute in complex. Loop in function runs from start index to i < end index. The end
 * index should therefore be for instance the size of the solute.
 * @param resolution Surface sites per atom. Default is 64.
 * @return Vector of surface sites of visible solute sites.
 */
std::vector<MolecularSurface::SurfaceSite> getVisibleMolecularSurface(const AtomCollection& molecule, int startIndex,
                                                                      int endIndex, int resolution = 64);
/**
 * @brief Helper function to write molecular surface into an xyz file with the surface points represented by H atoms.
 *
 * @param os The output stream to write to.
 * @param surface Vector of surface sites.
 */
void writeSurface(std::ostream& os, std::vector<MolecularSurface::SurfaceSite> surface);

} /* namespace MolecularSurface */
} /* namespace Utils */
} /* namespace Scine */

#endif // UTILSOS_MOLECULARSURFACE_H
