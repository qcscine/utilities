/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/GeometryOptimization/NtOptimizer.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Geometry/ElementInfo.h"
#include "Utils/GeometryOptimization/NtOptimizerSettings.h"
#include "Utils/Optimizer/LeastSquares/LevenbergMarquardt.h"
#include <math.h>
#include <Eigen/Dense>
#include <array>
#include <fstream>
#include <iostream>
#include <valarray>
namespace Scine {
namespace Utils {

int NtOptimizer::optimize(AtomCollection& atoms, Core::Log& log) {
  // Configure Calculator
  _calculator.setStructure(atoms);
  _calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);
  // Transform into internal coordinates
  Utils::PositionCollection coordinates = atoms.getPositions();
  const unsigned int nAtoms = atoms.size();
  int cycle = 0;
  _norms.clear();
  for (unsigned int loop = 0; loop < this->check.maxIter; loop++) {
    cycle++;
    // Calculate data
    Utils::Results results = _calculator.calculate("Geometry Optimization Cycle");
    if (!results.get<Property::SuccessfulCalculation>()) {
      throw std::runtime_error("Gradient calculation in geometry optimization failed.");
    }
    double value = results.get<Property::Energy>();
    auto gradients = results.get<Property::Gradients>();
    this->triggerObservers(cycle, value, Eigen::Map<const Eigen::VectorXd>(coordinates.data(), nAtoms * 3));
    // Evaluate additional force
    double rcNorm = this->updateGradients(atoms, value, gradients);
    _norms.push_back(rcNorm);
    _trajectory.push_back(coordinates);
    // Check convergence
    auto positions = atoms.getPositions();
    Eigen::Vector3d lhsCenter = Eigen::Vector3d::Zero();
    for (auto& l : this->lhsList) {
      lhsCenter += positions.row(l);
    }
    lhsCenter.array() /= this->lhsList.size();
    Eigen::Vector3d rhsCenter = Eigen::Vector3d::Zero();
    for (auto& r : this->rhsList) {
      rhsCenter += positions.row(r);
    }
    rhsCenter.array() /= this->rhsList.size();
    double c2c = (lhsCenter - rhsCenter).norm();
    if (this->attractive) {
      if (c2c < this->check.attractiveStop) {
        coordinates = this->extractTsGuess();
        atoms.setPositions(coordinates);
        _calculator.modifyPositions(coordinates);
        return cycle;
      }
      for (auto& l : this->lhsList) {
        for (auto& r : this->rhsList) {
          Eigen::Vector3d dir = positions.row(l) - positions.row(r);
          const double dist = dir.norm();
          const double r12cov =
              ElementInfo::covalentRadius(atoms.getElement(l)) + ElementInfo::covalentRadius(atoms.getElement(r));
          if (dist < this->check.attractiveStop * r12cov) {
            coordinates = this->extractTsGuess();
            atoms.setPositions(coordinates);
            _calculator.modifyPositions(coordinates);
            return cycle;
          }
        }
      }
    }
    else {
      bool done = true;
      for (auto& l : this->lhsList) {
        for (auto& r : this->rhsList) {
          Eigen::Vector3d dir = positions.row(l) - positions.row(r);
          const double dist = dir.norm();
          const double r12cov =
              ElementInfo::covalentRadius(atoms.getElement(l)) + ElementInfo::covalentRadius(atoms.getElement(r));
          if (dist < this->check.repulsiveStop * r12cov) {
            done = false;
            break;
          }
        }
      }
      if (c2c <= this->check.repulsiveStop) {
        done = false;
      }
      if (done) {
        coordinates = this->extractTsGuess();
        atoms.setPositions(coordinates);
        _calculator.modifyPositions(coordinates);
        return cycle;
      }
    }
    // Update positions (SD)
    if (this->transformCoordinates) {
      try {
        auto transformation = std::make_shared<InternalCoordinates>(atoms);
        auto internalCoordinates = transformation->coordinatesToInternal(coordinates);
        auto internalGradients = transformation->gradientsToInternal(gradients);
        internalCoordinates -= optimizer.factor * internalGradients;
        coordinates = transformation->coordinatesToCartesian(internalCoordinates);
      }
      catch (const InternalCoordinatesException& e) {
        coordinates -= optimizer.factor * gradients;
      }
    }
    else {
      coordinates -= optimizer.factor * gradients;
    }
    atoms.setPositions(coordinates);
    _calculator.modifyPositions(coordinates);
  }
  return cycle;
}

double NtOptimizer::updateGradients(const AtomCollection& atoms, const double& energy, GradientCollection& gradients) const {
  // Gather data
  auto positions = atoms.getPositions();
  const unsigned int nAtoms = atoms.size();
  const auto& lhs = this->lhsList;
  const auto& rhs = this->rhsList;

  if (lhs.size() == 0 || rhs.size() == 0)
    throw std::logic_error("Called NT optimization without atoms to apply a force to.");

  // Check lists
  for (auto& i : lhs) {
    if (i >= nAtoms) {
      throw std::logic_error("Index greater than number of atoms requested in NT atom list.");
    }
  }
  for (auto& i : rhs) {
    if (i >= nAtoms) {
      throw std::logic_error("Index greater than number of atoms requested in NT atom list.");
    }
  }

  // define reaction coordinate based in LHS and RHS list
  GradientCollection reactionCoordinate(gradients);
  reactionCoordinate.setZero();
  Eigen::Vector3d lhsCenter;
  Eigen::Vector3d rhsCenter;
  lhsCenter.setZero();
  rhsCenter.setZero();
  for (auto& l : lhs) {
    lhsCenter += positions.row(l);
  }
  lhsCenter.array() /= lhs.size();
  for (auto& r : rhs) {
    rhsCenter += positions.row(r);
  }
  rhsCenter.array() /= rhs.size();
  Eigen::VectorXd c2c = lhsCenter - rhsCenter;
  c2c /= c2c.norm();
  for (auto& l : lhs) {
    reactionCoordinate.row(l) = c2c;
  }
  for (auto& r : rhs) {
    reactionCoordinate.row(r) = -c2c;
  }
  // remove gradient along reaction coordinate
  double rcNorm = reactionCoordinate.norm();
  const double alignment = gradients.cwiseProduct(reactionCoordinate).sum() / rcNorm;
  GradientCollection gAlongRC = (alignment / rcNorm) * reactionCoordinate;
  //  GradientCollection gResidual = gradients set to zero for contrained atoms;
  GradientCollection gResidual = gradients;
  for (auto& l : lhs) {
    gResidual.row(l).setZero();
  }
  for (auto& r : rhs) {
    gResidual.row(r).setZero();
  }

  // replace gradient along reaction coordinate with fixed force
  if (this->attractive) {
    gradients = gResidual + 0.5 * totalForceNorm * reactionCoordinate;
  }
  else {
    gradients = gResidual - 0.5 * totalForceNorm * reactionCoordinate;
  }

  double sign = (alignment > 0.0) ? 1.0 : -1.0;

  //  std::ofstream outfile;
  //  outfile.open("nt.data.txt", std::ios_base::app);
  //  outfile << sign * gAlongRC.norm() << " " << energy << std::endl;

  return sign * gAlongRC.norm();
}

PositionCollection NtOptimizer::extractTsGuess() const {
  // Try to find a clear TS guess as indicated by a minus to plus zero pass
  //   Search from the back if the force is attractive
  if (this->attractive) {
    for (unsigned int i = _norms.size() - 2; i > 0; i--) {
      if (_norms[i] <= 0.0 && _norms[i + 1] > 0.0) {
        return _trajectory[(fabs(_norms[i]) < fabs(_norms[i + 1])) ? i : i + 1];
      }
    }
  }
  else {
    for (unsigned int i = 0; i < _norms.size() - 1; i++) {
      if (_norms[i] >= 0.0 && _norms[i + 1] < 0.0) {
        return _trajectory[(fabs(_norms[i]) < fabs(_norms[i + 1])) ? i : i + 1];
      }
    }
  }
  throw std::runtime_error("No transition state guess was found in Newton Trajectory scan.");
  return _trajectory[_trajectory.size() - 1];
};

void NtOptimizer::addSettingsDescriptors(UniversalSettings::DescriptorCollection& collection) const {
  // Not implemented yet
  throw std::runtime_error("You reached a function in the NTOptimizer that should not be called.");
};

void NtOptimizer::setSettings(const Settings& settings) {
  // optimizer.applySettings(settings);
  this->optimizer.factor = settings.getDouble(NtOptimizer::ntSdFactorKey);
  this->check.maxIter = settings.getInt(NtOptimizer::ntMaxIterKey);
  this->check.repulsiveStop = settings.getDouble(NtOptimizer::ntRepulsiveStopKey);
  this->check.attractiveStop = settings.getDouble(NtOptimizer::ntAttractiveStopKey);
  this->rhsList = settings.getIntList(NtOptimizer::ntRHSListKey);
  this->lhsList = settings.getIntList(NtOptimizer::ntLHSListKey);
  this->attractive = settings.getBool(NtOptimizer::ntAttractiveKey);
  this->totalForceNorm = settings.getDouble(NtOptimizer::ntTotalForceNormKey);
  this->transformCoordinates = settings.getBool(NtOptimizer::ntTransfromCoordinatesKey);
};

void NtOptimizer::applySettings(const Settings& settings) {
  this->setSettings(settings);
};

Settings NtOptimizer::getSettings() const {
  return NtOptimizerSettings(*this);
}

} // namespace Utils
} // namespace Scine
