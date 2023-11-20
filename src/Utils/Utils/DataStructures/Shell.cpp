/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Shell.h"
#include <utility>

namespace Scine {
namespace Utils {
namespace Integrals {

Shell::Shell() = default;

Shell::~Shell() = default;

Shell::Shell(std::vector<double> v_alpha, std::vector<double> v_coeffs, Displacement shift, std::size_t max_l, bool pure_solid)
  : _vecAlpha(std::move(v_alpha)),
    _vecCoeffs(std::move(v_coeffs)),
    _shift(std::move(shift)),
    _l(max_l),
    _pureSolid(pure_solid),
    _contrLegnth(_vecAlpha.size()) {
  if (_vecCoeffs.size() != _contrLegnth) {
    throw std::runtime_error("Length of coefficient vector is not eqial to length of alpha vector.");
  }
  _calcMaxLogCoeffs();
} // constructor

} // namespace Integrals
} // namespace Utils
} // namespace Scine
