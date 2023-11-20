/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_ELECTRONICOCCUPATIONGENERATOR_H
#define UTILS_ELECTRONICOCCUPATIONGENERATOR_H

namespace Scine {
namespace Utils {

class LcaoMethod;
namespace LcaoUtils {
class ElectronicOccupation;

/*!
 * This interface generates an ElectronicOccupation instance, from which the density matrix can be generated (in
 * combination with the molecular orbitals).
 */
class ElectronicOccupationGenerator {
 public:
  virtual ~ElectronicOccupationGenerator() = default;

  void setMethod(LcaoMethod* method);
  ElectronicOccupation generateOccupation();
  /*! Notify that a new SCF cycle has started, allows to reinitialize some members if necessary. */
  void newScfCycleStarted();

 protected:
  LcaoMethod* method_ = nullptr;

 private:
  virtual ElectronicOccupation generateOccupationImpl() = 0;
  virtual void newScfCycleStartedImpl() {
  }
  bool occupationFulfillsRequirements(const ElectronicOccupation& occupation);
};

} // namespace LcaoUtils

} // namespace Utils
} // namespace Scine
#endif // UTILS_ELECTRONICOCCUPATIONGENERATOR_H
