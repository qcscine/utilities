/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/GeometryOptimization/AfirOptimizerBase.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Geometry/ElementInfo.h"
#include <math.h>
#include <Eigen/Dense>
#include <array>

namespace Scine {
namespace Utils {

void AfirOptimizerBase::evaluateArtificialForces(const AtomCollection& atoms, double& energy,
                                                 GradientCollection& gradients) const {
  // clang-format off
  constexpr const std::array<double, 119> covalentRadii = {{ 0.0, // dummy
   32.0,                                                                                                                  46.0,
  133.0, 102.0,                                                                        85.0,  75.0,  71.0,  63.0,  64.0,  67.0,
  155.0, 139.0,                                                                       126.0, 116.0, 111.0, 103.0,  99.0,  96.0,
  196.0, 171.0, 148.0, 136.0, 134.0, 122.0, 119.0, 116.0, 111.0, 110.0, 112.0, 118.0, 124.0, 121.0, 121.0, 116.0, 114.0, 117.0,
  210.0, 185.0, 163.0, 154.0, 147.0, 138.0, 128.0, 125.0, 125.0, 120.0, 128.0, 136.0, 142.0, 140.0, 140.0, 136.0, 133.0, 131.0,
  232.0, 196.0,
                180.0, 163.0, 176.0, 174.0, 173.0, 172.0, 168.0, 169.0, 168.0, 167.0, 166.0, 165.0, 164.0, 170.0, 162.0,
                       152.0, 146.0, 137.0, 131.0, 129.0, 122.0, 123.0, 124.0, 133.0, 144.0, 144.0, 151.0, 145.0, 147.0, 142.0,
  223.0, 201.0,
                186.0, 175.0, 169.0, 170.0, 171.0, 172.0, 166.0, 166.0, 168.0, 168.0, 165.0, 167.0, 173.0, 176.0, 161.0,
                       157.0, 149.0, 143.0, 141.0, 134.0, 129.0, 128.0, 121.0, 122.0, 136.0, 143.0, 162.0, 175.0, 165.0, 157.0}};
  // clang-format on

  // Gather data
  auto positions = atoms.getPositions();
  auto elements = atoms.getElements();
  const unsigned int nAtoms = atoms.size();
  const auto& lhs = this->lhsList;
  const auto& rhs = this->rhsList;

  // Check lists
  for (auto& i : lhs) {
    if (i >= nAtoms) {
      throw std::logic_error("Index greater than number of atoms requested in AFIR atom list.");
    }
  }
  for (auto& i : rhs) {
    if (i >= nAtoms) {
      throw std::logic_error("Index greater than number of atoms requested in AFIR atom list.");
    }
  }

  // Clean initialization
  energy = 0.0;
  gradients.setZero();

  // Calculate some constants
  const double epsilon = 0.00038320319; // 1.0061 kJ/mol in a.u.
  const double r0 = 7.21195078;         // 3.8164 Angstrom in a.u.
  const double gamma = this->energyAllowance / 2625.5;
  const double tmp = 1.0 + sqrt(1.0 + gamma / epsilon);
  const double prefactor = (this->attractive ? 1.0 : -1.0) * gamma / ((pow(2.0, -1.0 / 6.0) - pow(tmp, -1.0 / 6.0)) * r0);
  const double gammaWeak = 10.0 / ((nAtoms - 1.0) * 2625.5);
  const double tmpWeak = 1.0 + sqrt(1.0 + gammaWeak / epsilon);
  const double prefactorWeak = gammaWeak / ((pow(2.0, -1.0 / 6.0) - pow(tmpWeak, -1.0 / 6.0)) * r0);

  if (this->lhsList.size() > 0 && this->rhsList.size() > 0) {
    // Generate distances and weights
    Eigen::MatrixXd distances = Eigen::MatrixXd::Zero(this->lhsList.size(), this->rhsList.size());
    Eigen::MatrixXd weights = Eigen::MatrixXd::Zero(this->lhsList.size(), this->rhsList.size());
    for (unsigned int i = 0; i < this->lhsList.size(); i++) {
      for (unsigned int j = 0; j < this->rhsList.size(); j++) {
        distances(i, j) = (positions.row(lhs[i]) - positions.row(rhs[j])).norm();
        // TODO add covalent radii
        double radi = covalentRadii[ElementInfo::Z(elements[lhs[i]])] * 0.0188971616463207;
        double radj = covalentRadii[ElementInfo::Z(elements[rhs[j]])] * 0.0188971616463207;
        weights(i, j) = (radi + radj) / distances(i, j);
      }
    }

    // Precalculate data and generate the energy contribution
    Eigen::MatrixXd wPow6 = Eigen::MatrixXd::Zero(this->lhsList.size(), this->rhsList.size());
    // The line below is not used due to some problems with it when linking
    //  against MKL
    // wPow6.array() = weights.array().pow(6);
    wPow6.array() = weights.array() * weights.array() * weights.array() * weights.array() * weights.array() * weights.array();
    const double pow6sum = wPow6.array().sum();
    const double pow5sum = (wPow6.array() * distances.array()).sum();
    energy += prefactor * (pow5sum / pow6sum);
    // Calculate the gradient contributions
    for (unsigned int i = 0; i < this->lhsList.size(); i++) {
      for (unsigned int j = 0; j < this->rhsList.size(); j++) {
        const double tmpPow6 = wPow6(i, j);
        const double tmpPow7 = wPow6(i, j) / distances(i, j);
        const double dEdr = ((6.0 * tmpPow7 * pow5sum) - (5.0 * tmpPow6 * pow6sum)) / (pow6sum * pow6sum);
        const auto v = positions.row(lhs[i]) - positions.row(rhs[j]);
        const double drdx = v[0] / distances(i, j);
        const double drdy = v[1] / distances(i, j);
        const double drdz = v[2] / distances(i, j);
        gradients(lhs[i], 0) += dEdr * drdx;
        gradients(lhs[i], 1) += dEdr * drdy;
        gradients(lhs[i], 2) += dEdr * drdz;
      }
    }
    for (unsigned int j = 0; j < this->rhsList.size(); j++) {
      for (unsigned int i = 0; i < this->lhsList.size(); i++) {
        const double tmpPow6 = wPow6(i, j);
        const double tmpPow7 = wPow6(i, j) / distances(i, j);
        const double dEdr = ((6.0 * tmpPow7 * pow5sum) - (5.0 * tmpPow6 * pow6sum)) / (pow6sum * pow6sum);
        const auto v = positions.row(rhs[j]) - positions.row(lhs[i]);
        const double drdx = v[0] / distances(i, j);
        const double drdy = v[1] / distances(i, j);
        const double drdz = v[2] / distances(i, j);
        gradients(rhs[j], 0) += dEdr * drdx;
        gradients(rhs[j], 1) += dEdr * drdy;
        gradients(rhs[j], 2) += dEdr * drdz;
      }
    }
    gradients *= prefactor;
  }
  /* This is what is denoted as '+' in AFIR1+ in:
   *  J. Chem. Theory Comput. ASAP (DOI: 10.1021/acs.jctc.8b01182)
   * https://pubs.acs.org/doi/10.1021/acs.jctc.8b01182
   */
  if (this->weak) {
    // Generate distances and weights
    Eigen::MatrixXd fullDistances = Eigen::MatrixXd::Zero(nAtoms, nAtoms);
    Eigen::MatrixXd fullWeights = Eigen::MatrixXd::Zero(nAtoms, nAtoms);
    for (unsigned int i = 0; i < nAtoms; i++) {
      for (unsigned int j = 0; j < i; j++) {
        fullDistances(i, j) = (positions.row(i) - positions.row(j)).norm();
        fullDistances(j, i) = fullDistances(i, j);
        double radi = covalentRadii[ElementInfo::Z(elements[i])] * 0.0188971616463207;
        double radj = covalentRadii[ElementInfo::Z(elements[j])] * 0.0188971616463207;
        fullWeights(i, j) = (radi + radj) / fullDistances(i, j);
        fullWeights(j, i) = fullWeights(i, j);
      }
    }

    // Precalculate data and generate the energy contribution
    Eigen::MatrixXd fullPow6 = Eigen::MatrixXd::Zero(nAtoms, nAtoms);
    // The line below is not used due to some problems with it when linking
    //  against MKL
    // fullPow6.array() = fullWeights.array().pow(6);
    fullPow6.array() = fullWeights.array() * fullWeights.array() * fullWeights.array() * fullWeights.array() *
                       fullWeights.array() * fullWeights.array();
    const double fullPow6sum = fullPow6.array().sum() / 2.0;
    const double fullPow5sum = (fullPow6.array() * fullDistances.array()).sum() / 2.0;
    energy += prefactorWeak * (fullPow5sum / fullPow6sum);
    // Calculate the gradient contributions
    auto weakGradients(gradients);
    weakGradients.setZero();
    for (unsigned int i = 0; i < nAtoms; i++) {
      for (unsigned int j = 0; j < nAtoms; j++) {
        if (i == j)
          continue;
        const double tmpPow6 = fullPow6(i, j);
        const double tmpPow7 = fullPow6(i, j) / fullDistances(i, j);
        const double dEdr = ((6.0 * tmpPow7 * fullPow5sum) - (5.0 * tmpPow6 * fullPow6sum)) / (fullPow6sum * fullPow6sum);
        const auto v = positions.row(i) - positions.row(j);
        const double drdx = v[0] / fullDistances(i, j);
        const double drdy = v[1] / fullDistances(i, j);
        const double drdz = v[2] / fullDistances(i, j);
        weakGradients(i, 0) += dEdr * drdx;
        weakGradients(i, 1) += dEdr * drdy;
        weakGradients(i, 2) += dEdr * drdz;
      }
    }
    gradients += prefactorWeak * weakGradients;
  }
}

} // namespace Utils
} // namespace Scine
