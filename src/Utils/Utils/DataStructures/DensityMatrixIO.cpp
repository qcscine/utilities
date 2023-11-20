/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DensityMatrixIO.h"
#include "DensityMatrix.h"
#include <fstream>

namespace Scine {
namespace Utils {

using namespace std;
void DensityMatrixIO::write(const std::string& filename, const DensityMatrix& m) {
  /*
   * Format : - restricted / unrestricted
   *          - number atomic orbitals
   *          - number alpha electrons
   *          - number beta electrons
   *          - one or two matrices, depending on restricted or unrestricted
   */

  ofstream fout;
  fout.open(filename, ios_base::out | ios_base::trunc | ios_base::binary);

  auto unrestricted = static_cast<int8_t>(m.unrestricted());
  fout.write(reinterpret_cast<char*>(&unrestricted), sizeof(int8_t)); // NOLINT

  auto nAOs = static_cast<int32_t>(m.restrictedMatrix().rows());
  fout.write(reinterpret_cast<char*>(&nAOs), sizeof(int32_t)); // NOLINT

  auto nAlpha = static_cast<int32_t>(m.numberElectronsInAlphaMatrix());
  auto nBeta = static_cast<int32_t>(m.numberElectronsInBetaMatrix());
  fout.write(reinterpret_cast<char*>(&nAlpha), sizeof(int32_t)); // NOLINT
  fout.write(reinterpret_cast<char*>(&nBeta), sizeof(int32_t));  // NOLINT

  if (unrestricted != 0) {
    const auto& alpha = m.alphaMatrix();
    const auto& beta = m.betaMatrix();
    fout.write(reinterpret_cast<const char*>(alpha.data()), nAOs * nAOs * sizeof(double)); // NOLINT
    fout.write(reinterpret_cast<const char*>(beta.data()), nAOs * nAOs * sizeof(double));  // NOLINT
  }
  else {
    const auto& matrix = m.restrictedMatrix();
    fout.write(reinterpret_cast<const char*>(matrix.data()), nAOs * nAOs * sizeof(double)); // NOLINT
  }
}

DensityMatrix DensityMatrixIO::read(const std::string& filename) {
  /*
   * Format : - restricted / unrestricted
   *          - number atomic orbitals
   *          - one or two matrices, depending on restricted or unrestricted
   */

  ifstream fin;
  fin.open(filename, ios_base::in | ios_base::binary);
  int8_t unrestricted8;
  int32_t nAOs, nAlpha, nBeta;
  fin.read(reinterpret_cast<char*>(&unrestricted8), sizeof(int8_t)); // NOLINT
  auto unrestricted = unrestricted8 != 0;
  fin.read(reinterpret_cast<char*>(&nAOs), sizeof(int32_t));   // NOLINT
  fin.read(reinterpret_cast<char*>(&nAlpha), sizeof(int32_t)); // NOLINT
  fin.read(reinterpret_cast<char*>(&nBeta), sizeof(int32_t));  // NOLINT

  DensityMatrix d;
  if (unrestricted) {
    Eigen::MatrixXd alpha(nAOs, nAOs), beta(nAOs, nAOs);
    fin.read(reinterpret_cast<char*>(alpha.data()), nAOs * nAOs * sizeof(double)); // NOLINT
    fin.read(reinterpret_cast<char*>(beta.data()), nAOs * nAOs * sizeof(double));  // NOLINT
    d.setDensity(std::move(alpha), std::move(beta), nAlpha, nBeta);
  }
  else {
    Eigen::MatrixXd matrix(nAOs, nAOs);
    fin.read(reinterpret_cast<char*>(matrix.data()), nAOs * nAOs * sizeof(double)); // NOLINT
    d.setDensity(std::move(matrix), nAlpha + nBeta);
  }

  return d;
}

} // namespace Utils
} // namespace Scine
