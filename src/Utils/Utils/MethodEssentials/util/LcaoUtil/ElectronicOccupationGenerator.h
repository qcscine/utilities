/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_ELECTRONICOCCUPATIONGENERATOR_H
#define UTILS_ELECTRONICOCCUPATIONGENERATOR_H

namespace Scine {
namespace Utils {

class LCAOMethod;
namespace LcaoUtil {
class ElectronicOccupation;

/*!
 * This interface generates an ElectronicOccupation instance, from which the density matrix can be generated (in
 * combination with the molecular orbitals).
 */
class ElectronicOccupationGenerator {
 public:
  virtual ~ElectronicOccupationGenerator() = default;

  void setMethod(LCAOMethod* method);
  ElectronicOccupation generateOccupation();
  /*! Notify that a new SCF cycle has started, allows to reinitialize some members if necessary. */
  void newScfCycleStarted();

 protected:
  LCAOMethod* method_ = nullptr;

 private:
  virtual ElectronicOccupation generateOccupationImpl() = 0;
  virtual void newScfCycleStartedImpl() {
  }
  bool occupationFulfillsRequirements(const ElectronicOccupation& occupation);
};

} // namespace LcaoUtil

} // namespace Utils
} // namespace Scine
#endif // UTILS_ELECTRONICOCCUPATIONGENERATOR_H
