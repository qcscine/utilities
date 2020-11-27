/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_DENSITYMATRIXIO_H
#define UTILS_DENSITYMATRIXIO_H

#include <string>

namespace Scine {
namespace Utils {

class DensityMatrix;

/*!
 * Class to write density matrices to disk or read them from disk, in binary format.
 * TODO: save space and only write and read triangular matrix (since symmetric)?
 * TODO: Make this class employ EigenMatrixIO? NB: then the memory layout would change
 */
class DensityMatrixIO {
 public:
  /*! Write the density matrix to disk. */
  static void write(const std::string& filename, const DensityMatrix& m);

  /*! Read a density matrix from the disk. */
  static DensityMatrix read(const std::string& filename);
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_DENSITYMATRIXIO_H
