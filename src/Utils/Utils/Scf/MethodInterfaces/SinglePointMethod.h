/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_SINGLEPOINTMETHOD_H
#define UTILS_SINGLEPOINTMETHOD_H

#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/Math/AtomicSecondDerivativeCollection.h>
#include <Utils/Math/DerivOrderEnum.h>
#include <Utils/Math/FullSecondDerivativeCollection.h>
#include <Utils/Typenames.h>
#include <vector>

namespace Scine {
namespace Core {
struct Log;
} // namespace Core

namespace Utils {

class AtomCollection;

/*!
 * This class is the base class of any method to be used in
 * Real-time quantum chemistry.
 * It contain declarations for variables methods concerning energy and gradients.
 */

class SinglePointMethod {
 public:
  explicit SinglePointMethod(Utils::DerivativeOrder maximalOrder);
  virtual ~SinglePointMethod() = default;

  /*! Performs one single-point calculation of the energy.
   * \param derivativesOrder which derivative to calculate up to, if possible. */
  virtual void calculate(Utils::Derivative d, Core::Log& log) = 0;

  void setAtomCollection(const Utils::AtomCollection& structure);
  /*! Set up the structure from the positions and element types.
      \param positions Vector of positions, in bohr. */
  void initializeStructure(const Utils::ElementTypeCollection& elements, const Utils::PositionCollection& positions);
  /*! Set up the structure from the element types. */
  void initializeStructure(const Utils::ElementTypeCollection& elements);
  /*! Update the positions. Do not use to initialize the structure!. */
  void setPositions(Utils::PositionCollection positions);
  /*! Get the positions (bohr). */
  const Utils::PositionCollection& getPositions() const;
  /*! Get the gradients (Hartree / bohr). */
  const Utils::GradientCollection& getGradients() const;
  /*! Get the second derivatives with respect to the single atomic nuclei (Hartree / bohr^2). */
  const Utils::AtomicSecondDerivativeCollection& getAtomicSecondDerivatives() const;
  /*! Get the full second derivatives containers. */
  const Utils::FullSecondDerivativeCollection& getFullSecondDerivatives() const;
  /*! Get the elements. */
  const Utils::ElementTypeCollection& getElementTypes() const;

  /*! Returns true if the method is able to calculate second derivatives. */
  bool canCalculateSecondDerivatives() const;

  /*! Returns the lower-triangular bond-order matrix. */
  const Utils::BondOrderCollection& getBondOrderCollection() const;
  void setBondOrderCollection(const Utils::BondOrderCollection& B);

  int getNumberAtoms() const;

  /*! Function returning the energy in Hartree. */
  double getEnergy() const;
  /*! Set the energy in Hartree. */
  void setEnergy(double energy);

  const std::vector<double>& getAtomicCharges() const;
  double getAtomicCharge(int index) const;
  void setAtomicCharges(const std::vector<double>& charges);

 protected:
  /*! Resets the size of gradients_, positions_, bondOrderMatrix_ and atomicCharges_. To be called during initialization
   * of the specific method. */
  void resizeRealTimeMethodMembers();

  Utils::BondOrderCollection bondOrders_;
  Utils::ElementTypeCollection elementTypes_;
  Utils::PositionCollection positions_;
  const Utils::DerivativeOrder maximalCalculableDerivativeOrder_;
  Utils::GradientCollection gradients_;
  Utils::AtomicSecondDerivativeCollection secondDerivatives_;
  Utils::FullSecondDerivativeCollection fullSecondDerivatives_;
  std::vector<double> atomicCharges_;
  double energy_;
};

inline int SinglePointMethod::getNumberAtoms() const {
  return static_cast<int>(elementTypes_.size());
}

} // namespace Utils
} // namespace Scine
#endif // UTILS_SINGLEPOINTMETHOD_H
