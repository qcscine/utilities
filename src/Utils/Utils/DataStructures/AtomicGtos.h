/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_ATOMICGTOS_H
#define UTILS_ATOMICGTOS_H

#include "Utils/DataStructures/Gtf.h"
#include "Utils/DataStructures/GtoExpansion.h"
#include <boost/optional.hpp>
#include <boost/variant.hpp>
#include <map>
#include <string>

namespace Scine {
namespace Utils {

//! Class containing GTO expansions for one atom.
struct AtomicGtos {
  boost::optional<GtoExpansion> s;
  boost::optional<GtoExpansion> p;
  boost::optional<GtoExpansion> d;

  std::map<std::string, std::vector<Gtf>> getGtfs() {
    std::map<std::string, std::vector<Gtf>> out;
    if (s) {
      out.insert(std::make_pair("s", s->gtfs));
    }
    if (p) {
      out.insert(std::make_pair("p", p->gtfs));
    }
    if (d) {
      out.insert(std::make_pair("d", d->gtfs));
    }
    return out;
  };
  /* variant type is necessary, because the first entry is the angular momentum and all other entries are tuples
   * describing the Gtfs */
  using nwChemFormatAngularMomentum = std::vector<boost::variant<int, std::pair<double, double>>>;
  using nwChemFormat = std::vector<nwChemFormatAngularMomentum>;
  nwChemFormat getNWChemFormat() {
    nwChemFormat out;
    nwChemFormatAngularMomentum values;
    if (s) {
      values.push_back(0); // angular momentum
      for (const auto& gtf : s->gtfs) {
        double exp = gtf.exponent;
        double coeff = gtf.coefficient;
        values.push_back(std::make_pair(exp, coeff));
      }
      out.push_back(values);
    }
    values.clear();
    if (p) {
      values.push_back(1); // angular momentum
      for (const auto& gtf : p->gtfs) {
        double exp = gtf.exponent;
        double coeff = gtf.coefficient;
        values.push_back(std::make_pair(exp, coeff));
      }
      out.push_back(values);
    }
    values.clear();
    if (d) {
      values.push_back(2); // angular momentum
      for (const auto& gtf : d->gtfs) {
        double exp = gtf.exponent;
        double coeff = gtf.coefficient;
        values.push_back(std::make_pair(exp, coeff));
      }
      out.push_back(values);
    }
    return out;
  };
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_ATOMICGTOS_H
