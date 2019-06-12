/**
 * @file
 * @brief A file containing definitions of the results
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_RESULTS_H
#define UTILS_RESULTS_H

#include "Utils/MethodEssentials/util/DipoleMatrix.h"
#include "Utils/MethodEssentials/util/SpinAdaptedMatrix.h"
#include "Utils/Typenames.h"
#include <exception>
#include <memory>
#include <string>

namespace Scine {
namespace Utils {

// Forward declarations
class BondOrderCollection;

/**
 * Exception thrown if a property is requested and it is not stored.
 */
class PropertyNotPresentException : public std::exception {
 public:
  const char* what() const noexcept final {
    return "Property desired not present in results.";
  }
};

/**
 * Class for the properties obtained in a single-point calculation.
 * To obtain the properties:
 * - "get" methods return references to the results (also allows simple copy)
 * - "take" methods move the results (which will not be present in Results afterwards).
 *
 * Implementation note:
 * - New voices can be added to the Result. Minimally a set, a get and a has method must be present.
 * - Use pImpl idiom to hide boost::optional dependency
 */
class Results {
 public:
  Results();
  ~Results();
  Results(Results&& rhs) noexcept;
  Results& operator=(Results&& rhs) noexcept;
  Results(const Results& rhs);
  Results& operator=(const Results& rhs);

  bool hasDescription() const;
  void setDescription(std::string description);
  const std::string& getDescription() const;

  bool hasEnergy() const;
  void setEnergy(double e);
  double getEnergy() const;

  bool hasGradients() const;
  void setGradients(GradientCollection gradients);
  const GradientCollection& getGradients() const;
  GradientCollection takeGradients();

  bool hasHessian() const;
  void setHessian(HessianMatrix hessian);
  const HessianMatrix& getHessian() const;
  HessianMatrix takeHessian();

  bool hasDipole() const;
  void setDipole(Dipole dipole);
  const Dipole& getDipole() const;

  bool hasDipoleGradient() const;
  void setDipoleGradient(DipoleGradient dipoleGradient);
  const DipoleGradient& getDipoleGradient() const;
  DipoleGradient takeDipoleGradient();

  bool hasAODipoleMatrix() const;
  void setAODipoleMatrix(DipoleMatrix dipole);
  const DipoleMatrix& getAODipoleMatrix() const;
  DipoleMatrix takeAODipoleMatrix();

  bool hasMODipoleMatrix() const;
  void setMODipoleMatrix(DipoleMatrix dipole);
  const DipoleMatrix& getMODipoleMatrix() const;
  DipoleMatrix takeMODipoleMatrix();

  bool hasOneElectronMatrix() const;
  void setOneElectronMatrix(Eigen::MatrixXd oneElectronMatrix);
  const Eigen::MatrixXd& getOneElectronMatrix() const;
  Eigen::MatrixXd takeOneElectronMatrix();

  bool hasTwoElectronMatrix() const;
  void setTwoElectronMatrix(SpinAdaptedMatrix oneElectronMatrix);
  const SpinAdaptedMatrix& getTwoElectronMatrix() const;
  SpinAdaptedMatrix takeTwoElectronMatrix();

  bool hasBondOrders() const;
  void setBondOrders(BondOrderCollection);
  const BondOrderCollection& getBondOrders() const;
  BondOrderCollection takeBondOrders();

 private:
  struct Impl; // pImpl idiom
  std::unique_ptr<Impl> pImpl_;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_RESULTS_H
