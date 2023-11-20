/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILSOS_SOLUTESOLVENTCOMPLEX_H
#define UTILSOS_SOLUTESOLVENTCOMPLEX_H

#include "Utils/Geometry.h"
#include "Utils/Solvation/MolecularSurface.h"
#include "Utils/Typenames.h"

namespace Scine {

namespace Core {
struct Log;
}

namespace Utils {
/**
 * @brief A tool for systematically solvating one solvate
 *
 * Namespace to set up solute solvent clusters as starting structures for calculations
 */
namespace SoluteSolventComplex {
/**
 * @brief Structure of settings for solvent placement
 *
 * @param resolution Number of surface sites per atom
 * @param solventOffset Initial distance to arrange solvent molecule.
 * @param maxDistance Maximum distance below it should be attempted to add a solvent molecule before the next surface
 * point is tested out.
 * @param stepSize Increment of solvent offset if no solvent molecule could be added.
 * @param numRotamers Number of rotations to be tried for adding solvent.
 * @param strategicSolv Bool to turn on or off the solvationStrategy function. Reduces the number of surface
 * recalculations.
 * @param coverageThreshold Ratio of visible solvent sites to be covered by solvents per total number of visible solvent
 * sites for the given solute. Once the threshold is reached the shell is considered complete.
 */
struct SolventPlacementSettings {
  int resolution = 18;
  double solventOffset = 0.0;
  double maxDistance = 5.0;
  double stepSize = 1.0;
  int numRotamers = 3;
  bool strategicSolv = true;
  double coverageThreshold = 0.85;
};
/**
 * @brief Add systematically a number of one type of solvent to solute
 *
 * @param soluteComplex AtomCollection of the solute (solute complex, if already partially solvated)
 * @param soluteSize Size of the solute, assuming the solute starts with index 0.
 * @param solvent AtomCollection of the solvent.
 * @param numSolvents Number of solvent molecules which should be added.
 * @param seed Seed for shuffling the surface sites and choosing the surface site of the solvent.
 * @param placementSettings Settings for the placement of one solvent.
 *
 * @return A solvent shell vector. Each entry is one vector of atom collections, containing the solvent molecules of
 * this shell.
 */
std::vector<std::vector<AtomCollection>> solvate(const AtomCollection& soluteComplex, int soluteSize,
                                                 const AtomCollection& solvent, int numSolvents, int seed,
                                                 SolventPlacementSettings placementSettings);
/**
 * @brief Add systematically a number of different types of solvent with the given ratio to the solute.
 *
 * @param soluteComplex AtomCollection of the solute (solute complex, if already partially solvated)
 * @param soluteSize Size of the solute, assuming the solute starts with index 0.
 * @param solvents A vector of AtomCollections of different solvents.
 * @param solventRatios A vector of ratios of the different solvents.
 * @param numSolvents Number of solvent molecules which should be added.
 * @param seed Seed for shuffling the surface sites, shuffling the occurrence of the different solvents and choosing
 * the surface site of the solvent.
 * @param placementSettings Settings for the placement of one solvent.
 *
 * @return A tuple of a solvent shell vector and a solvent shell indices vector (to determine where which solvent is).
 * Each entry of the solvent shell vector is one vector of atom collections, containing the solvent molecules of
 * this shell.
 * Each entry of the solvent shell indices vector is one vector of indices, containing the solvent indices of the
 * solvents given in the original input.
 */
std::tuple<std::vector<std::vector<AtomCollection>>, std::vector<std::vector<int>>>
solvateMix(const AtomCollection& soluteComplex, int soluteSize, const std::vector<AtomCollection>& solvents,
           const std::vector<int>& solventRatios, int numSolvents, int seed, SolventPlacementSettings placementSettings);
/**
 * @brief Add number of shells of one type of solvent to solute.
 *
 * @param soluteComplex AtomCollection of the solute (solute complex, if already partially solvated)
 * @param soluteSize Size of the solute, assuming the solute starts with index 0.
 * @param solvent AtomCollection of solvent.
 * @param numShells Number of shells to be completed by one type of solvent. One shell means that no surface point
 * is 'visible' anymore.
 * @param seed Seed for shuffling the surface sites and choosing the surface site of the solvent
 * @param placementSettings Settings for the placement of one solvent.
 *
 * @return A solvent shell vector. Each entry is one vector of atom collections, containing the solvent molecules of
 * this shell.
 */
std::vector<std::vector<AtomCollection>> solvateShells(const AtomCollection& soluteComplex, int soluteSize,
                                                       const AtomCollection& solvent, int numShells, int seed,
                                                       SolventPlacementSettings placementSettings);
/**
 * @brief Add number of shells of different solvents to solute.
 *
 * @param soluteComplex AtomCollection of the solute (solute complex, if already partially solvated)
 * @param soluteSize Size of the solute, assuming the solute starts with index 0.
 * @param solvents A vector of AtomCollections of different solvents.
 * @param solventRatios A vector of ratios of the different solvents.
 * @param numShells Number of shells to be completed by solvents with given ratio. One shell means that no surface point
 * is 'visible' anymore.
 * @param seed Seed for shuffling the surface sites, shuffling the occurrence of the different solvents and choosing
 * the surface site of the solvent.
 * @param placementSettings Settings for the placement of one solvent.
 *
 * @return A tuple of a solvent shell vector and a solvent shell indices vector (to determine where which solvent is).
 * Each entry of the solvent shell vector is one vector of atom collections, containing the solvent molecules of
 * this shell.
 * Each entry of the solvent shell indices vector is one vector of indices, containing the solvent indices of the
 * solvents given in the original input.
 */
std::tuple<std::vector<std::vector<AtomCollection>>, std::vector<std::vector<int>>>
solvateShellsMix(const AtomCollection& soluteComplex, int soluteSize, const std::vector<AtomCollection>& solvents,
                 const std::vector<int>& solventRatios, int numShells, int seed, SolventPlacementSettings placementSettings);
/**
 * @brief Add systematically a number of solvents and a number of solvent shells to solute
 *
 * @param soluteComplex AtomCollection of the solute (solute complex, if already partially solvated)
 * @param soluteSize Size of the solute, assuming the solute starts with index 0.
 * @param solvent AtomCollection of solvent.
 * @param numShells Number of shells to be completed. One shell means that no surface point is 'visible' anymore.
 * @param numSolvents Number of solvent molecules which should be added.
 * @param seed Seed for shuffling the surface sites and choosing the surface site of the solvent
 * @param placementSettings Settings for the placement of one solvent.
 *
 * @return A tuple of a solvent shell vector and a solvent shell indices vector (to determine where which solvent is).
 * Each entry of the solvent shell vector is one vector of atom collections, containing the solvent molecules of
 * this shell.
 * Each entry of the solvent shell indices vector is one vector of indices, containing the solvent indices of the
 * solvents given in the original input.
 */
std::tuple<std::vector<std::vector<AtomCollection>>, std::vector<std::vector<int>>>
solvate(const AtomCollection& soluteComplex, int soluteSize, const std::vector<AtomCollection>& solvents,
        const std::vector<int>& solventRatios, int numSolvents, int numShells, int seed,
        SolventPlacementSettings placementSettings);
/**
 * @brief Create a vector of indices of the solvent molecules.
 *
 * @param numSolvents Number of solvent molecules which should be added.
 * @param ratios A vector of ratios of the different solvents.
 * @param differentSolvents Number of different solvents.
 *
 * @return A vector with indices of solvents.
 */
std::vector<int> getSolventIndices(int numSolvents, const std::vector<int>& ratios, unsigned long int differentSolvents);

/**
 * @brief Analyze given solute-solvent complex to return solvent shell vector. Basically reverse solvate function,
 * giving the solvent shell vector of a given cluster.
 *
 * @param complex Solute solvent cluster to be analyzed.
 * @param soluteSize Original size of the solute.
 * @param solventSizeVector Vector containing the solvent sizes in order.
 * @param resolution Number of surface sites per atom.
 * @param strategicSolv Bool to turn on or off the solvationStrategy function. Reduces the number of surface
 * recalculations when a solvent has been added.
 * @param coverageThreshold Ratio of visible solvent sites to be covered by solvents per total number of visible solvent
 * sites for the given solute. Once the threshold is reached the shell is considered complete.

 * @return A solvent shell vector in a solute solvent complex. Each entry is one vector of atom collections, containing
 * the solvent molecules of this shell. Threshold allows of more loose definition of solvent shell vector.
 */
std::vector<std::vector<AtomCollection>>
giveSolventShellVector(const AtomCollection& complex, int soluteSize, const std::vector<int>& solventSizeVector,
                       int resolution, Core::Log& log, bool strategicSolv = true, double coverageThreshold = 1.0);
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
 * @brief Flatten the solvent shell indices vector to one dimension
 *
 * @param solventShellIndices A vector of vectors of indices.
 *
 * @return A vector of indices to track which input solvent is where.
 */
std::vector<int> mergeSolventShellIndices(const std::vector<std::vector<int>>& solventShellIndices);

/**
 * @brief Transform solvent shell vector to solvent size vector. Needed for giveSolventShellVector calculation.
 *
 * @param shellVector Solvent shell vector containing a vector of shells containing a vector of solvent molecules
 * @return A one dimensional vector containing the size of all solvents in order of their addition by the solvate
 * function.
 */
std::vector<int> transferSolventShellVector(const std::vector<std::vector<AtomCollection>>& shellVector);
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
PositionCollection arrange(const Position& surfPoint1, const Position& surfNormal1, const Position& surfPoint2,
                           const Position& surfNormal2, const PositionCollection& atoms2Position, double distance);
/**
 * @brief Add additive to given complex at given surface site of the complex.
 *
 * @param complex AtomCollection of complex to which additive should be added. Called by reference.
 * @param additive AtomCollection of additive supposed to be added to complex.
 * @param complexSurfSite SurfaceSite of complex along which additive should be arranged.
 * @param additiveSurfSite SurfaceSite of additive along which additive should be arranged.
 * @param minDistance Minimal distance between two surface sites attempted to add.
 * @param maxDistance Maximum distance between two surface sites attempted to add.
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
