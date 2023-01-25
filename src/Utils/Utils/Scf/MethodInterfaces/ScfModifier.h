/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_SCFMODIFIER_H
#define UTILS_SCFMODIFIER_H

namespace Scine {
namespace Utils {

class ScfMethod;

/**
 * @class ScfModifier ScfModifier.h
 *
 * @brief Base class for SCF modifiers.
 *
 * SCF modifiers allow to act in between the different steps of a SCF calculations.
 * The functions it defines are anchors and are called after the corresponding steps of a SCF calculation.
 * \sa ScfMethod
 */

class ScfModifier {
 public:
  virtual ~ScfModifier() = default;

  /**
   * @brief Set the SCF method; the access points of the SCF modifier will be called from this SCF method.
   *
   * This function can be overridden to perform further checks / ...
   */
  virtual void setMethod(ScfMethod* method) {
    m = method;
  }

  /**
   * @brief Initialize the SCF Modifier; to be overridden if necessary by the derived class.
   */
  virtual void initialize() {
  }

  /**
   * @brief Gives information about compatibility between ScfMethod and ScfModifier; to be overridden if necessary by
   * derived classes.
   * @return true if the ScfMethod and ScfModifier are compatible, true by default.
   */
  virtual bool isValid() {
    return true;
  }

  /**
   * @brief Function called after the calculation of the overlap matrix.
   */
  virtual void onOverlapCalculated() {
  }

  /**
   * @brief Function called at the start of every SCF iteration.
   */
  virtual void onIterationStart() {
  }

  /**
   * @brief Function called after the calculation of the Fock matrix.
   */
  virtual void onFockCalculated() {
  }

  /**
   * @brief Function called after the SCF eigenvalue problem has been solved.
   */
  virtual void onGEPSolved() {
  }

  /**
   * @brief Function called after the density matrix has been calculated.
   */
  virtual void onDensityCalculated() {
  }

  /**
   * @brief Function called at the very end of a SCF calculation.
   */
  virtual void onCalculationFinalized() {
  }

 protected:
  ScfMethod* m = nullptr;
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_SCFMODIFIER_H
