/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_PARTIALHESSIAN_H
#define UTILS_PARTIALHESSIAN_H

#include "Utils/Geometry/AtomCollection.h"
#include <Eigen/Core>
#include <utility>
#include <vector>

namespace Scine {
namespace Utils {

/*!
 * Class defining a partial Hessian for use in embedding calculations.
 * It is simply a container for the matrix and the indices of the atoms in the supersystem.
 */
class PartialHessian {
 public:
  using Matrix = Eigen::MatrixXd;

  PartialHessian(Matrix hessian, std::vector<int> indices)
    : _hessian(std::move(hessian)), _indices(std::move(indices)) {
    if (_hessian.rows() != _hessian.cols()) {
      throw std::runtime_error("The given Hessian is not square.");
    }
    if (_hessian.rows() != static_cast<long>(3 * _indices.size())) {
      throw std::runtime_error("The given Hessian does not match the given indices.");
    }
  }

  PartialHessian(const PartialHessian& other) = default;
  PartialHessian& operator=(const PartialHessian& other) = default;

  inline const Matrix& getMatrix() const {
    return _hessian;
  }

  inline const std::vector<int>& getIndices() const {
    return _indices;
  }

  inline int getNumberOfAtoms() const {
    return static_cast<int>(_indices.size());
  }

  inline AtomCollection getPartialAtoms(const ElementTypeCollection& elements,
                                        const PositionCollection& positionCollection) const {
    auto atoms = AtomCollection(elements, positionCollection);
    return getPartialAtoms(atoms);
  }

  inline AtomCollection getPartialAtoms(const AtomCollection& atoms) const {
    int nAtoms = atoms.size();
    AtomCollection partialAtoms;
    for (int i : _indices) {
      if (nAtoms < i) {
        throw std::runtime_error("The given atoms cannot be the correct super system, there are too few atoms.");
      }
      partialAtoms.push_back(atoms[i]);
    }
    return partialAtoms;
  }

 private:
  Matrix _hessian;
  std::vector<int> _indices;
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_PARTIALHESSIAN_H
