/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_STRUCTUREDEPENDENTINITIALIZER_H
#define UTILS_STRUCTUREDEPENDENTINITIALIZER_H

#include <Utils/Typenames.h>
#include <vector>

namespace Scine {
namespace Utils {

class AtomsOrbitalsIndexes;

/*!
 * Interface for method-specific settings which must be updated f.i. when the molecular structure changes.
 */
class StructureDependentInitializer {
 public:
  virtual ~StructureDependentInitializer() = default;

  /*! Initialize the method <b>after</b> the parameters have been set or loaded. */
  virtual void initialize(const Utils::ElementTypeCollection& elements) = 0;

  virtual AtomsOrbitalsIndexes getAtomsOrbitalsIndexes() const = 0;
  virtual unsigned getNumberElectronsForUnchargedSpecies() const = 0;
  virtual std::vector<double> getCoreCharges() const = 0;
  virtual bool unrestrictedCalculationPossible() const = 0;
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_STRUCTUREDEPENDENTINITIALIZER_H
