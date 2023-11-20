/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_GTOEXPANSION_H
#define UTILS_GTOEXPANSION_H

#include <Utils/DataStructures/Gtf.h>
#include <vector>

namespace Scine {
namespace Utils {

/*!
 * The class GTOExpansion is the container for the coefficients of a
 * STO-nG expansion.
 */
struct GtoExpansion {
  int angularMomentum = 0;
  std::vector<Gtf> gtfs;

  int nAOs() const {
    return 2 * angularMomentum + 1;
  };
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_GTOEXPANSION_H
