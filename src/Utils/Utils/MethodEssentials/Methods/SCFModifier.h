/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_SCFMODIFIER_H
#define UTILS_SCFMODIFIER_H

namespace Scine {
namespace Utils {

class SCFMethod;

/*!
 * Base class for SCF modifiers.
 * SCF modifiers allow to act in between the different steps of a SCF calculations.
 * The functions it defines are anchors and are called after the corresponding steps of a SCF calculation.
 * \sa SCFMethod
 */

class SCFModifier {
 public:
  virtual ~SCFModifier() = default;

  /*! Set the method its access points will be called from.
      This function can be overriden to perform further checks / ... */
  virtual void setMethod(SCFMethod* method) {
    m = method;
  }

  /*! Initialize the SCF Modifier, to be overridden if necessary by the derived class. */
  virtual void initialize() {
  }

  /*! Gives information about compatibility between SCFMethod and SCFModifier.
      To be overriden if necessary by derived classes.
      \return true if the SCFMethod and SCFModifier are compatible, true by default. */
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
  SCFMethod* m = nullptr;
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_SCFMODIFIER_H
