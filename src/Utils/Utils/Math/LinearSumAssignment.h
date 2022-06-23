/**
 * @file LinearSumAssignment.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_LINEAR_SUM_ASSIGNMENT_H
#define UTILS_LINEAR_SUM_ASSIGNMENT_H

#include <Eigen/Core>
#include <exception>
#include <map>
#include <numeric>
#include <vector>

namespace Scine {
namespace Utils {
namespace LinearSumAssignment {

enum class Type { Minimize, Maximize };

struct DualVariables {
  std::vector<double> u;
  std::vector<double> v;
};

constexpr const double infinity = std::numeric_limits<double>::max();
constexpr const int notAssigned = -1;

struct InvalidCostMatrixDimensionsException : public std::exception {
  auto what() const noexcept -> const char* final;
};
struct InvalidCostMatrixException : public std::exception {
  auto what() const noexcept -> const char* final;
};

template<Type type>
void prepareCostMatrix(Eigen::MatrixXd& /*costMatrix*/) {
}
template<>
inline void prepareCostMatrix<Type::Maximize>(Eigen::MatrixXd& costMatrix) {
  costMatrix.array() *= -1;
  costMatrix.array() -= costMatrix.minCoeff();
}

/**
 * @brief Class handling the augmentation of the minimal-cost path in Linear Sum Assignment.
 * @class PathAugmenter @file LinearSumAssignment.h
 * Finding the shortest augmenting path is done according to
 *
 * David F. Crouse, On Implementing 2D Rectangular Assignment Algorithms,
 * IEEE Transactions on Aerospace and Electronic Systems, 2016, 52, 1679.
 * doi: 10.1109/TAES.2016.140952
 *
 * This class is instantiated once every iteration of the Linear Sum Assignment solver.
 * In total, it is instantiated a number of times equal to the minimum of the number of
 * rows and cols in the cost matrix.
 * Obviously, some vectors could be initialized outside this class, but honestly,
 * this will never have a real impact on the runtime and I like it better like this
 * because it is way more readable.
 *
 * This clas, after the path augmentation, has the responsibility to update the
 * values of the dual variables with the knowledge of the minimal cost path.
 * The current sink of the path can also be queried.
 */
struct PathAugmenter {
  PathAugmenter(const Eigen::MatrixXd& costMatrix);

  /**
   * @brief Augments the current minimal cost path.
   */
  void augment(int currentRow, const std::vector<int>& row4col, std::vector<int>& path, const DualVariables& duals);

  /**
   * @brief Updates the dual variables U and V.
   */
  auto updateDualVariables(int currentRow, const std::vector<int>& col4row, DualVariables currentDuals) -> DualVariables;

  /**
   * @brief Queries what column is the last of the current path.
   */
  auto sink() const -> int;

 private:
  const Eigen::MatrixXd& costMatrix_;
  std::vector<double> shortestPathCosts_;
  std::vector<bool> SR_;
  std::vector<bool> SC_;
  int sink_{-1};
  double minVal_{0};
};

/**
 * @brief Class doing a Linear Sum Assignment.
 * @class Solver @file LinearSumAssignment.h
 *
 * The Linear Sum Assignment algorithm is implemented according to
 *
 * David F. Crouse, On Implementing 2D Rectangular Assignment Algorithms,
 * IEEE Transactions on Aerospace and Electronic Systems, 2016, 52, 1679.
 * doi: 10.1109/TAES.2016.140952
 *
 * This solver handles quadratic and rectangular cost matrices, and can both
 * solve minimization and maximization problems.
 *
 * @tparam type A LinearSumAssignment::Type, flags whether minimation or maximization is
 *              required
 *
 * Usage:
 * @code{.cpp}
 * Eigen::MatrixXd costMatrix = Eigen::MatrixXd::Identity(20, 30);
 * LinearSumAssignment::Solver<LinearSumAssignment::Type::Maximize> solver(costMatrix);
 *
 * // Obtain a map with as key a row index and as value the assigned column index
 * auto result = solver();
 * // Obtain the assignment cost:
 * double cost = solver.cost();
 * @endcode
 *
 */
template<Type type = Type::Minimize>
struct Solver {
  /**
   * @brief Constructor. Takes the cost matrix, checks its sanity and transposes if nCols < nRows.
   *
   *  For maximizations, the cost matrix C is modified as:
   *  # cost -> score
   *  C *= -1;
   *  # make matrix positive
   *  C -= min(C)
   */
  Solver(Eigen::MatrixXd costMatrix);
  /**
   * @brief Solve the linear sum assignment problem given a cost matrix.
   */
  auto operator()() -> std::map<int, int>;
  /**
   * @brief Returns the cost (or score for maximizations) of the assignment.
   */
  auto cost() const -> double;
  /**
   * @brief Returns the result map.
   */
  auto result() const -> const std::map<int, int>&;

 private:
  void prepareInput();
  void augmentSolution(int currentRow, int sink);
  Eigen::MatrixXd costMatrix_;
  bool transposed_{false};
  std::vector<int> col4row_;
  std::vector<int> row4col_;
  std::vector<int> path_;
  DualVariables dualVariables_;
  std::map<int, int> result_;
};

template<Type type>
Solver<type>::Solver(Eigen::MatrixXd costMatrix)
  : costMatrix_(std::move(costMatrix)), transposed_(costMatrix_.rows() > costMatrix_.cols()) {
  if (costMatrix_.rows() == 0 || costMatrix_.cols() == 0) {
    throw InvalidCostMatrixDimensionsException();
  }
  if (transposed_) {
    costMatrix_.transposeInPlace();
  }
  prepareCostMatrix<type>(costMatrix_);
}

template<Type type>
void Solver<type>::prepareInput() {
  col4row_ = std::vector<int>(costMatrix_.rows(), notAssigned);
  row4col_ = std::vector<int>(costMatrix_.cols(), notAssigned);
  path_ = std::vector<int>(costMatrix_.cols(), notAssigned);
  dualVariables_ = {std::vector<double>(costMatrix_.rows(), 0), std::vector<double>(costMatrix_.cols(), 0)};
}

template<Type type>
auto Solver<type>::operator()() -> std::map<int, int> {
  prepareInput();

  for (int currentRow = 0; currentRow < costMatrix_.rows(); ++currentRow) {
    PathAugmenter augmenter(costMatrix_);

    augmenter.augment(currentRow, row4col_, path_, dualVariables_);
    dualVariables_ = augmenter.updateDualVariables(currentRow, col4row_, dualVariables_);
    augmentSolution(currentRow, augmenter.sink());
  }

  result_ = {};
  for (int i = 0; i < costMatrix_.rows(); ++i) {
    result_.insert(transposed_ ? std::pair<int, int>{col4row_[i], i} : std::pair<int, int>{i, col4row_[i]});
  }

  return result_;
}

template<Type type>
void Solver<type>::augmentSolution(int currentRow, int sink) {
  int i = notAssigned;
  int j = sink;
  do {
    i = path_[j];
    row4col_[j] = i;
    std::swap(col4row_[i], j);

  } while (i != currentRow);
}

template<Type type>
auto Solver<type>::cost() const -> double {
  return std::accumulate(result_.begin(), result_.end(), 0.0, [&](double current, const std::pair<int, int>& element) -> double {
    return current + costMatrix_(element.first, element.second);
  });
}

template<Type type>
auto Solver<type>::result() const -> const std::map<int, int>& {
  return result_;
}
} // namespace LinearSumAssignment
} // namespace Utils
} // namespace Scine

#endif
