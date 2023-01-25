/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "ShellPairs.h"

namespace Scine {
namespace Utils {
namespace Integrals {

ShellPair::ShellPair(const Shell& shell, bool calculateCauchySchwarzFactor)
  : calculateCauchySchwarzFactor_(calculateCauchySchwarzFactor), shell_(std::make_unique<Shell>(shell)) {
}

ShellPair::iterator ShellPair::begin() {
  return shellPairDataVector_.begin();
}

ShellPair::const_iterator ShellPair::begin() const {
  return shellPairDataVector_.begin();
}

ShellPair::iterator ShellPair::end() {
  return shellPairDataVector_.end();
}

ShellPair::const_iterator ShellPair::end() const {
  return shellPairDataVector_.end();
}

ShellPairData& ShellPair::front() {
  return shellPairDataVector_.front();
}

const ShellPairData& ShellPair::front() const {
  return shellPairDataVector_.front();
}

int ShellPair::size() const {
  return shellPairDataVector_.size();
}

ShellPairData& ShellPair::operator[](int index) {
  return shellPairDataVector_[index];
}

const ShellPairData& ShellPair::operator[](int index) const {
  return shellPairDataVector_[index];
}

ShellPairData& ShellPair::at(int index) {
  return shellPairDataVector_.at(index);
}

const ShellPairData& ShellPair::at(int index) const {
  return shellPairDataVector_.at(index);
}

void ShellPairs::resize(int newSize) {
  shellPairs_.resize(newSize);
}

ShellPairs::iterator ShellPairs::begin() {
  return shellPairs_.begin();
}
ShellPairs::const_iterator ShellPairs::begin() const {
  return shellPairs_.begin();
}
ShellPairs::iterator ShellPairs::end() {
  return shellPairs_.end();
}
ShellPairs::const_iterator ShellPairs::end() const {
  return shellPairs_.end();
}
int ShellPairs::size() const {
  return shellPairs_.size();
}
bool ShellPairs::empty() const {
  return shellPairs_.empty();
}
ShellPair& ShellPairs::operator[](int index) {
  return shellPairs_[index];
}
const ShellPair& ShellPairs::operator[](int index) const {
  return shellPairs_[index];
}
ShellPair& ShellPairs::at(int index) {
  return shellPairs_.at(index);
}
const ShellPair& ShellPairs::at(int index) const {
  return shellPairs_.at(index);
}

const Shell& ShellPair::getShell() const {
  return *shell_;
}

ShellPairs::~ShellPairs() = default;

} // namespace Integrals
} // namespace Utils
} // namespace Scine
