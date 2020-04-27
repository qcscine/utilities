/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_LEVENBERGMARQUARDT_H_
#define UTILS_LEVENBERGMARQUARDT_H_

#include "UpdateFunctionManagerBase.h"
#include <Utils/Optimizer/Optimizer.h>

namespace Scine {
namespace Utils {
/**
 * @class LevenbergMarquardt LevenbergMarquardt.h
 * @brief Levenberg-Marquardt (LM) algorithm for a least squares optimization.
 */
class LevenbergMarquardt : public Optimizer {
 public:
  static constexpr const char* calculateCovarianceMatrixKey = "calculate_covariance_matrix";
  static constexpr const char* maxFuncEvalKey = "max_function_evaluations";
  /// @brief Default constructor.
  LevenbergMarquardt() = default;
  /**
   * @brief This function performs the optimization.
   * @param parameters The parameters to be optimized in a least squares way.
   * @param updateFunctionManager The update function manager provided to calculate errors and the Jacobian.
   */
  void optimize(Eigen::VectorXd& parameters, UpdateFunctionManagerBase& updateFunctionManager);
  /**
   * @struct Functor required to be specified for the Eigen LM algorithm.
   */
  struct LMFunctor {
    /**
     * @brief Constructor from the update function manager.
     */
    explicit LMFunctor(UpdateFunctionManagerBase& updateFunctionManager);
    /**
     * @brief Calculates the errors for a given set of parameters.
     */
    int operator()(const Eigen::VectorXd& parameters, Eigen::VectorXd& fvec) const;
    /**
     * @brief Calculates the Jacobian of the errors for a given set of parameters.
     */
    int df(const Eigen::VectorXd& parameters, Eigen::MatrixXd& fjac) const;
    /**
     * @brief Number of data points, i.e. values.
     */
    int m;
    /**
     * @brief The number of parameters, i.e. inputs.
     */
    int n;
    /**
     * @brief Getter for values.
     * \return Number of values m.
     */
    int values() const;
    /**
     * @brief Getter for inputs
     * @return Number of inputs n.
     */
    int inputs() const;
    /**
     * @brief The underlying update function manager.
     */
    UpdateFunctionManagerBase& updateFunctionManager_;
  };
  /**
   * @brief Adds all relevant options to the given UniversalSettings::DescriptorCollection
   *        thus expanding it to include the least squares options.
   * @param collection The DescriptorCollection to which new fields shall be added.
   */
  void addSettingsDescriptors(UniversalSettings::DescriptorCollection& collection) const final {
    UniversalSettings::BoolDescriptor calculateCovarianceMatrixDescriptor(
        "Calculate a covariance matrix for the parameters along with the optimization.");
    calculateCovarianceMatrixDescriptor.setDefaultValue(true);
    collection.push_back(calculateCovarianceMatrixKey, std::move(calculateCovarianceMatrixDescriptor));

    UniversalSettings::IntDescriptor maxFuncEvalDescriptor(
        "Sets the maximum number of function evaluations during an optimization.");
    maxFuncEvalDescriptor.setDefaultValue(0); // Because it will then not be passed on to Eigen, which will then use its
                                              // default value.
    collection.push_back(maxFuncEvalKey, std::move(maxFuncEvalDescriptor));
  };
  /**
   * @brief Updates the LM's options with those values given in the Settings.
   * @param settings The settings to update the option of the LM algorithm with.
   */
  void applySettings(const Settings& settings) final {
    calculateCovarianceMatrix = settings.getBool(calculateCovarianceMatrixKey);
    maxFuncEval = settings.getInt(maxFuncEvalKey);
  };

  /**
   * @brief Maximum number of function evaluations during an optimization
   *
   *        The default value is set to zero here, but this value is only
   *        passed on to the internal Eigen parameter if it set to be greater
   *        than zero!
   */
  int maxFuncEval = 0;
  /// @brief Decides whether a covariance matrix should be calculated for the parameters.
  bool calculateCovarianceMatrix = true;
  /// @brief Getter for the covariance matrix.
  const Eigen::MatrixXd& getCovarianceMatrix();

 private:
  Eigen::MatrixXd covarianceMatrix_;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_LEVENBERGMARQUARDT_H_
