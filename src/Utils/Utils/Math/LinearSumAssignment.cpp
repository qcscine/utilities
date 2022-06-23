/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "LinearSumAssignment.h"

namespace Scine {
namespace Utils {
namespace LinearSumAssignment {
auto InvalidCostMatrixDimensionsException::what() const noexcept -> const char* {
  return "0 rows or columns in input cost matrix in Linear Sum Assignment";
}
auto InvalidCostMatrixException::what() const noexcept -> const char* {
  return "Only infinite costs detected in cost matrix in Linear Sum Assignment.";
}

PathAugmenter::PathAugmenter(const Eigen::MatrixXd& costMatrix) : costMatrix_(costMatrix) {
  shortestPathCosts_ = std::vector<double>(costMatrix.cols(), infinity);
  SR_ = std::vector<bool>(costMatrix.rows(), false);
  SC_ = std::vector<bool>(costMatrix.cols(), false);
}

void PathAugmenter::augment(int currentRow, const std::vector<int>& row4col, std::vector<int>& path, const DualVariables& duals) {
  assert(path.size() == row4col.size() && path.size() == SC_.size());
  int i = currentRow;
  sink_ = notAssigned;
  while (sink_ == notAssigned) {
    SR_[i] = true;

    for (int j = 0; j < int(SC_.size()); ++j) {
      double val = minVal_ + costMatrix_(i, j) - duals.u[i] - duals.v[j];
      if (!SC_[j] && shortestPathCosts_[j] > val) {
        path[j] = i;
        shortestPathCosts_[j] = val;
      }
    }
    int j = notAssigned;
    minVal_ = infinity;
    for (int h = 0; h < int(SC_.size()); ++h) {
      if (!SC_[h] && (shortestPathCosts_[h] < minVal_ || ((shortestPathCosts_[h] == minVal_ && row4col[h] == notAssigned)))) {
        j = h;
        minVal_ = shortestPathCosts_[h];
      }
    }
    if (minVal_ == infinity) {
      throw InvalidCostMatrixException();
    }
    SC_[j] = true;

    if (row4col[j] == notAssigned) {
      sink_ = j;
    }
    else {
      i = row4col[j];
    }
  }
}

auto PathAugmenter::updateDualVariables(int currentRow, const std::vector<int>& col4row, DualVariables currentDuals)
    -> DualVariables {
  currentDuals.u[currentRow] += minVal_;
  for (int i = 0; i < int(SR_.size()); ++i) {
    if (SR_[i] && i != currentRow) {
      currentDuals.u[i] += minVal_ - shortestPathCosts_[col4row[i]];
    }
  }
  for (int j = 0; j < int(SC_.size()); ++j) {
    if (SC_[j]) {
      currentDuals.v[j] += shortestPathCosts_[j] - minVal_;
    }
  }
  return currentDuals;
}

auto PathAugmenter::sink() const -> int {
  return sink_;
}
} // namespace LinearSumAssignment
} // namespace Utils
} // namespace Scine
