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

/*!
 * Base class for SCF modifiers.
 * SCF modifiers allow to act in between the different steps of a SCF calculations.
 * The functions it defines are anchors and are called after the corresponding steps of a SCF calculation.
 * \sa ScfMethod
 */

class ScfModifier {
 public:
  virtual ~ScfModifier() = default;

  /*! Set the method its access points will be called from.
      This function can be overriden to perform further checks / ... */
  virtual void setMethod(ScfMethod* method) {
    m = method;
  }

  /*! Initialize the SCF Modifier, to be overridden if necessary by the derived class. */
  virtual void initialize() {
  }

  /*! Gives information about compatibility between ScfMethod and ScfModifier.
      To be overriden if necessary by derived classes.
      \return true if the ScfMethod and ScfModifier are compatible, true by default. */
  virtual bool isValid() {
    return true;
  }

  // access points
  virtual void onOverlapCalculated() {
  }
  virtual void onIterationStart() {
  }
  virtual void onFockCalculated() {
  }
  virtual void onGEPSolved() {
  }
  virtual void onDensityCalculated() {
  }
  virtual void onCalculationFinalized() {
  }

 protected:
  ScfMethod* m = nullptr;
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_SCFMODIFIER_H
