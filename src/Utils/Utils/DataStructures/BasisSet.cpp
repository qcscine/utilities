/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "BasisSet.h"

namespace Scine {
namespace Utils {
namespace Integrals {

int BasisSet::_idProvider = 0;

std::vector<long> BasisSet::shellToAtom(const AtomCollection& atoms) const {
  return shellToAtom(*this, atoms, false);
}

std::vector<std::vector<long>> BasisSet::atomToShell(const AtomCollection& atoms) const {
  return atomToShell(atoms, *this);
}

std::vector<long> BasisSet::shellToAtom(const BasisSet& basis, const AtomCollection& atoms, bool throw_if_no_match = false) {
  std::vector<long> result;
  result.reserve(basis.size());
  for (const auto& shell : basis) {
    auto a = std::find_if(atoms.begin(), atoms.end(), [&shell](const Atom& a) {
      return shell.getShift()[0] == a.getPosition()[0] && shell.getShift()[1] == a.getPosition()[1] &&
             shell.getShift()[2] == a.getPosition()[2];
    });
    const auto found_match = (a != atoms.end());
    if (throw_if_no_match && !found_match) {
      throw std::logic_error("shell2atom: no matching atom found");
    }
    result.push_back(found_match ? a.get() - atoms.begin().get() : -1);
  }
  return result;
}

std::vector<std::vector<long>> BasisSet::atomToShell(const AtomCollection& atoms, const BasisSet& basis) {
  std::vector<std::vector<long>> result;
  result.resize(atoms.size());
  size_t iatom = 0;
  for (const auto& a : atoms) {
    auto s = basis.begin();
    while (s != basis.end()) {
      s = std::find_if(s, basis.end(), [&a](const Shell& s) {
        return s.getShift()[0] == a.getPosition()[0] && s.getShift()[1] == a.getPosition()[1] &&
               s.getShift()[2] == a.getPosition()[2];
      });
      if (s != basis.end()) {
        result[iatom].push_back(s - basis.begin());
        ++s;
      }
    }
    ++iatom;
  }
  return result;
}

const std::shared_ptr<ShellPairs>& BasisSet::getShellPairs() const {
  return _shellPairs;
}

void BasisSet::setShellPairs(std::shared_ptr<ShellPairs>& shellPairs) {
  _shellPairsEvaluated = true;
  _shellPairs = std::move(shellPairs);
}

void BasisSet::deleteShellPairs() {
  _shellPairsEvaluated = false;
  _shellPairs.reset();
}

bool BasisSet::areShellPairsEvaluated() const {
  return _shellPairsEvaluated;
}

auto BasisSet::append(const BasisSet& rhs) -> void {
  this->insert(this->end(), rhs.begin(), rhs.end());

  for (auto const& atom : rhs.getAtoms()) {
    _atoms.push_back(atom);
  }
}

const AtomCollection& BasisSet::getAtoms() const {
  return _atoms;
}

auto BasisSet::uncontract() -> void {
  BasisSet backUp = *this;

  (*this).resize(0);

  for (const auto& shell : backUp) {
    for (const auto& alpha : shell.getVecAlpha()) {
      bool alreadyInBasis = false;
      for (auto& shellAlreadyInBasis : (*this)) {
        if (std::abs((shell.getShift() - shellAlreadyInBasis.getShift()).norm()) < 1e-10 &&
            std::abs(alpha - shellAlreadyInBasis.getVecAlpha().at(0)) < 1e-10) {
          alreadyInBasis = true;
        }
      }
      if (!alreadyInBasis) {
        (*this).push_back(Shell({alpha}, {1.0}, shell.getShift(), shell.l(), shell.isPureSolid()));
      }
    }
  }
}

} // namespace Integrals
} // namespace Utils
} // namespace Scine
