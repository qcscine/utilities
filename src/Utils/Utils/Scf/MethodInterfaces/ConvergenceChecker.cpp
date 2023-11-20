/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ConvergenceChecker.h"
#include "ScfMethod.h"

namespace Scine {
namespace Utils {

ConvergenceChecker::ConvergenceChecker() {
  Thresholds t{1e-7, 1e-5};
  set(t);
}

ConvergenceChecker::ConvergenceChecker(Thresholds thresholds) {
  set(thresholds);
}

ConvergenceChecker::ConvergenceChecker(const ConvergenceChecker& rhs) {
  set(rhs.thresholds_);
  this->converged_ = rhs.converged_;
}
ConvergenceChecker& ConvergenceChecker::operator=(const ConvergenceChecker& rhs) {
  set(rhs.thresholds_);
  this->converged_ = rhs.converged_;
  return *this;
}

void ConvergenceChecker::set(ConvergenceChecker::Thresholds thresholds) {
  ConvergenceCheckerContainer types;
  thresholds_ = std::move(thresholds);
  if (thresholds_.energyDifference) {
    types.insert({ConvergenceType::EnergyDifference,
                  std::make_unique<ScfEnergyConvergenceChecker>(*thresholds_.energyDifference)});
  }
  if (thresholds_.densityRmsd) {
    types.insert({ConvergenceType::DensityRMSD, std::make_unique<ScfDensityConvergenceChecker>(*thresholds_.densityRmsd)});
  }
  convergenceTypes_ = std::move(types);
}

auto ConvergenceChecker::get() const -> Thresholds {
  return thresholds_;
}

auto ConvergenceChecker::getCurrentValues() const -> std::vector<boost::optional<double>> {
  std::vector<boost::optional<double>> values;
  values.reserve(convergenceTypes_.size());
  for (const auto& checker : convergenceTypes_) {
    values.push_back(checker.second->current());
  }
  return values;
}

void ConvergenceChecker::update(const ScfMethod& m) {
  for (auto& checker : convergenceTypes_) {
    checker.second->update(m);
  }
}

auto ConvergenceChecker::converged() const -> bool {
  for (const auto& checker : convergenceTypes_) {
    if (!checker.second->check()) {
      return false;
    }
  }
  // if no checking is given, then always return false:
  // scf goes on until maxiterations
  return !convergenceTypes_.empty();
}

auto ConvergenceChecker::getNames() const -> std::vector<std::string> {
  std::vector<std::string> names;
  names.reserve(convergenceTypes_.size());
  for (const auto& checker : convergenceTypes_) {
    names.push_back(checker.second->name());
  }
  return names;
}

auto ScfEnergyConvergenceChecker::check() const -> bool {
  // oldEnergy_ is set and threshold is respected
  return oldEnergy_ && *current_ <= threshold_;
}

void ScfEnergyConvergenceChecker::update(const ScfMethod& m) {
  oldEnergy_ = newEnergy_;
  newEnergy_ = m.getElectronicEnergy();
  if (oldEnergy_) {
    current_ = std::abs(*oldEnergy_ - *newEnergy_);
  }
}

auto ScfDensityConvergenceChecker::check() const -> bool {
  // size equal and threshold respected
  return (newMatrix_.size() == oldMatrix_.size()) && *current_ <= threshold_;
}

void ScfDensityConvergenceChecker::update(const ScfMethod& m) {
  std::swap(oldMatrix_, newMatrix_);
  newMatrix_ = m.getDensityMatrix().restrictedMatrix();
  if (newMatrix_.size() == oldMatrix_.size()) {
    current_ = (newMatrix_ - oldMatrix_).norm();
  }
}

} // namespace Utils
} // namespace Scine
