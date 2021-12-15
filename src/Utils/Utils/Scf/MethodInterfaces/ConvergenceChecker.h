/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_CONVERGENCECHECKER_H
#define UTILS_CONVERGENCECHECKER_H

#include <Eigen/Core>
#include <boost/optional.hpp>
#include <map>
#include <memory>
#include <vector>

namespace Scine {
namespace Utils {

class ScfMethod;

enum class ConvergenceType : unsigned { EnergyDifference, DensityRMSD };

struct ScfConvergenceChecker {
  ScfConvergenceChecker() = default;
  ScfConvergenceChecker(double threshold) : threshold_(threshold) {
  }
  virtual ~ScfConvergenceChecker() = default;

  /*! Interface method for the convergence checking */
  virtual auto check() const -> bool = 0;
  /*! Interface method for the updating the quantity in the convergence checking */
  virtual void update(const ScfMethod& /*m*/) {
  }
  /*! Name of the covergence criterion, important for printing. Less than 18 chars */
  virtual auto name() const -> std::string = 0;
  /*! Set the threshold under which the convergence criterion is satisfied.
      Corresponds to the rmsd in the density matrix change. */
  void setThreshold(double v) {
    threshold_ = v;
  }
  /*! Get the currently applied threshold. */
  auto getThreshold() const -> double {
    return threshold_;
  }

  /*! Get the current criterion value, for insance, the current energy error. */
  auto current() const -> boost::optional<double> {
    return current_;
  }

 protected:
  double threshold_;
  boost::optional<double> current_;
};

struct ScfEnergyConvergenceChecker final : public ScfConvergenceChecker {
 public:
  constexpr static const ConvergenceType Type = ConvergenceType::EnergyDifference;

  ScfEnergyConvergenceChecker(double threshold) : ScfConvergenceChecker(threshold) {
  }
  ~ScfEnergyConvergenceChecker() final = default;
  auto check() const -> bool final;
  void update(const ScfMethod& m) final;
  auto name() const -> std::string final {
    return "|Energy Difference|";
  }

 private:
  boost::optional<double> oldEnergy_;
  boost::optional<double> newEnergy_;
};

struct ScfDensityConvergenceChecker final : public ScfConvergenceChecker {
 public:
  constexpr static const ConvergenceType Type = ConvergenceType::DensityRMSD;
  ScfDensityConvergenceChecker(double threshold) : ScfConvergenceChecker(threshold) {
  }
  ~ScfDensityConvergenceChecker() final = default;
  auto check() const -> bool final;
  void update(const ScfMethod& m) final;
  auto name() const -> std::string final {
    return "Density RMSD";
  }

 private:
  Eigen::MatrixXd oldMatrix_, newMatrix_;
};

/**
 * @brief A Class handling the checking of scf convergence.
 * @class ConvergenceChecker @file ConvergenceChecker.h
 * Right now, energy convergence and density rmsd are implemented.
 *
 * Depending on what Thresholds are activated in set(...), the
 * right convergence criteria are checked
 * The default constructor sets the thresholds of the energy difference to
 * 1e-7 and the density rmsd to 1e-5.
 */
class ConvergenceChecker {
 public:
  using ConvergenceCheckerContainer = std::map<ConvergenceType, std::unique_ptr<ScfConvergenceChecker>>;

  struct Thresholds {
    boost::optional<double> energyDifference;
    boost::optional<double> densityRmsd;
  };

  ConvergenceChecker();
  explicit ConvergenceChecker(Thresholds thresholds);
  ~ConvergenceChecker() = default;
  ConvergenceChecker(const ConvergenceChecker& rhs);
  ConvergenceChecker(ConvergenceChecker&& rhs) noexcept = default;
  ConvergenceChecker& operator=(const ConvergenceChecker& rhs);
  ConvergenceChecker& operator=(ConvergenceChecker&& rhs) noexcept = default;

  /**
   * Set the threshold under which the convergence criterion is satisfied.
   * Corresponds to the energy difference between 2 iterations and/or the
   * rmsd in the density matrix change.
   * If no threshold is set, all checks return false (never converge).
   */
  void set(Thresholds thresholds);
  /**
   * Get the currently applied threshold.
   */
  auto get() const -> Thresholds;
  /**
   * Get the current values.
   */
  auto getCurrentValues() const -> std::vector<boost::optional<double>>;
  /**
   * Gets the current density matrix from the assigned SCF method.
   *  Performs no check of convergence and is useful if the starting density matrix in the SCF cycle is not the last one
   * of the previous cycle.
   * Having an update function avoids needing ScfMethod in the argument (non copyable)
   */
  void update(const ScfMethod& m);

  /**
   * Check if the convergence criterion is met.
   *   Implicitly calls the underlying convergece checkers.
   */
  auto converged() const -> bool;

  /**
   * Return the names of the convergence checkers in ConvergenceCheckerContainer
   * This is useful for printing.
   */
  auto getNames() const -> std::vector<std::string>;

 private:
  Thresholds thresholds_;
  ConvergenceCheckerContainer convergenceTypes_;
  bool converged_{false};
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_CONVERGENCECHECKER_H
