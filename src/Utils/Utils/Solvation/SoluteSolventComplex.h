/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILSOS_SOLUTESOLVENTCOMPLEX_H
#define UTILSOS_SOLUTESOLVENTCOMPLEX_H

#include "Utils/Geometry.h"
#include "Utils/Solvation/MolecularSurface.h"
#include "Utils/Typenames.h"

namespace Scine {
namespace Utils {
/**
 * @brief A tool for systematically solvating one solvate
 *
 * Namespace to set up solute solvent clusters as starting structures for calculations
 */
namespace SoluteSolventComplex {
/**
 * @brief Add systematically a number of solvents to solute
 *
 * @param soluteComplex AtomCollection of the solute (solute complex, if already partially solvated)
 * @param soluteSize Size of the solute, assuming the solute starts with index 0.
 * @param solvent AtomCollection of solvent.
 * @param numSolvents Number of solvent molecules which should be added.
 * @param seed Seed for shuffeling the surface sites and choosing the surface site of the solvent
 * @param resolution Number of surface sites per atom
 * @param solventOffset Initial distance to arrange solvent molecule.
 * @param maxDistance Maximum distance below it should be attempted to add a solvent molecule before the next surface
 * point is tested out.
 * @param stepSize Increament of solvent offset if no solvent molcule could be added.
 * @param numRotamers Number of rotations to be tried for adding solvent.
 * @param strategicSolv Bool to turn on or off the solvationStrategy function. Reduces the number of surface
 * recalculations.
 * @return A solvent shell vector. Each entry is one vector of atom collections, containg the solvent molecules of this
 * shell.
 */
std::vector<std::vector<AtomCollection>> solvate(const AtomCollection& soluteComplex, int soluteSize,
                                                 const AtomCollection& solvent, int numSolvents, int seed,
                                                 int resolution = 32, double solventOffset = 0.0, double maxDistance = 10.0,
                                                 double stepSize = 0.25, int numRotamers = 3, bool strategicSolv = false);
/**
 * @brief Add number of solvent shells to solute.
 *
 * @param soluteComplex AtomCollection of the solute (solute complex, if already partially solvated)
 * @param soluteSize Size of the solute, assuming the solute starts with index 0.
 * @param solvent AtomCollection of solvent.
 * @param numShells Number of shells to be completed. One shell means that no surface point is 'visible' anymore.
 * @param seed Seed for shuffeling the surface sites and choosing the surface site of the solvent
 * @param resolution Number of surface sites per atom
 * @param solventOffset Initial distance to arrange solvent molecule.
 * @param maxDistance Maximum distance below it should be attempted to add a solvent molecule before the next surface
 * point is tested out.
 * @param stepSize Increament of solvent offset if no solvent molcule could be added.
 * @param numRotamers Number of rotations to be tried for adding solvent.
 * @param strategicSolv Bool to turn on or off the solvationStrategy function. Reduces the number of surface
 * recalculations.
 * @return A solvent shell vector. Each entry is one vector of atom collections, containg the solvent molecules of this
 * shell.
 */
std::vector<std::vector<AtomCollection>>
solvateShells(const AtomCollection& soluteComplex, int soluteSize, const AtomCollection& solvent, int numShells,
              int seed, int resolution = 32, double solventOffset = 0.0, double maxDistance = 10.0,
              double stepSize = 0.25, int numRotamers = 3, bool strategicSolv = false);
/**
 * @brief Add systematically a number of solvents and a number of solvent shells to solute
 *
 * @param soluteComplex AtomCollection of the solute (solute complex, if already partially solvated)
 * @param soluteSize Size of the solute, assuming the solute starts with index 0.
 * @param solvent AtomCollection of solvent.
 * @param numShells Number of shells to be completed. One shell means that no surface point is 'visible' anymore.
 * @param numSolvents Number of solvent molecules which should be added.
 * @param seed Seed for shuffeling the surface sites and choosing the surface site of the solvent
 * @param resolution Number of surface sites per atom
 * @param solventOffset Initial distance to arrange solvent molecule.
 * @param maxDistance Maximum distance below it should be attempted to add a solvent molecule before the next surface
 * point is tested out.
 * @param stepSize Increament of solvent offset if no solvent molcule could be added.
 * @param numRotamers Number of rotations to be tried for adding solvent.
 * @param strategicSolv Bool to turn on or off the solvationStrategy function. Reduces the number of surface
 * recalculations.
 * @return A solvent shell vector. Each entry is one vector of atom collections, containg the solvent molecules of this
 * shell.
 */
std::vector<std::vector<AtomCollection>> solvate(const AtomCollection& soluteComplex, int soluteSize,
                                                 const AtomCollection& solvent, int numSolvents, int numShells,
                                                 int seed, int resolution, double solventOffset, double maxDistance,
                                                 double stepSize, int numRotamers, bool strategicSolv);
/**
 * @brief Merge a vector of atom collections to one atom collection.
 *
 * @param atomCollList A vector of AtomCollections
 * @return One AtomCollection, ordered as in the given vector.
 */
AtomCollection mergeAtomCollectionVector(const std::vector<AtomCollection>& atomCollList);
/**
 * @brief Merge a vector of a vector of atom collections (solvent shell vector) to one atom collection.
 *
 * @param atomCollList A vector of a vector of AtomCollections.
 * @return One AtomCollection, ordered as in the given vector.
 */
AtomCollection mergeSolventShellVector(const std::vector<std::vector<AtomCollection>>& shellVector);
/**
 * @brief Solvation strategy for faster building of solute - solvent complexes.
 * @param numberSurfPoints Number of visible surface points of the solute.
 * @return Number of solvents which should be added without recalculating the visible surface.
 */
int solvationStrategy(int numberSurfPoints);
/**
 * @brief Check if two atom collections overlap with their VdW spheres.
 * @param atoms1 First atom collection of interest.
 * @param atoms2 Second atom collection of interest.
 * @return Returns bool. True if the two collections do not overlap, False if they do overlap.
 */
bool checkDistances(const AtomCollection& atoms1, const AtomCollection& atoms2);
/**
 * @brief Arrange one atom collection such that the two positions given face each other.
 *
 * First position is of static atom collection and second atom collection is arranged relative to
 * it in such a way that the given surface points are set in the given distance and are facing each other.
 *
 * @param surfPoint1 Position of point at first atom collection.
 * @param surfNormal1 Position of normalized normal vector of surfPoint1.
 * @param surfPoint2 Position of point at second atom collection.
 * @param surfNormal2 Position of normalized normal vector of surfPoint2.
 * @param atoms2Position Position collection of atoms which should be arranged.
 * @param distance Targeted distance between surfPoint1 and surfPoint2.
 * @return Position collection of arranged atoms2.
 */
PositionCollection arrange(Position surfPoint1, Position surfNormal1, Position surfPoint2, Position surfNormal2,
                           const PositionCollection& atoms2Position, double distance);
/**
 * @brief Add additive to given complex at given surface site of the complex.
 *
 * @param complex AtomCollection of complex to which additive should be added. Called by reference.
 * @param additive AtomCollection of additive supposed to be added to complex.
 * @param complexSurfSite SurfaceSite of complex along which additive should be arranged.
 * @param additiveSurfSite SurfaceSite of additive along which additive should be arranged.
 * @param minDistance Minimal distance between two surface sites attempted to add.
 * @param maxDistance Maxium distance between two surface sites attempted to add.
 * @param incrementDistance Distance increment, if addition was unsuccessful and distance is still smaller than
 * maxDistance.
 * @param numRotationAttempts Number of rotations of additive to be tried to add additive.
 * @return Bool, true if addition was successful, false if addition was not successful; Complex gets manipulated as well.
 */
bool add(AtomCollection& complex, const AtomCollection& additive, const MolecularSurface::SurfaceSite& complexSurfSite,
         const MolecularSurface::SurfaceSite& additiveSurfSite, double minDistance, double maxDistance,
         double incrementDistance = 0.25, int numRotationAttempts = 3);

} /* namespace SoluteSolventComplex */

} /* namespace Utils */
} /* namespace Scine */

#endif // UTILSOS_SOLUTESOLVENTCOMPLEX_H
