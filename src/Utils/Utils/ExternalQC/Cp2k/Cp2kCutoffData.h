/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_EXTERNALQC_CP2KCUTOFFDATA_H
#define UTILS_EXTERNALQC_CP2KCUTOFFDATA_H

#include <cmath>
#include <utility>
#include <vector>

namespace Scine {
namespace Utils {

class AtomCollection;
class PropertyList;
class Settings;

namespace ExternalQC {

/**
 * @class Cp2kGridData Cp2kGridData.h
 * @brief This class holds data about the multigrid composition for the two possible cutoffs in CP2K
 */
class Cp2kCutoffData {
 public:
  Cp2kCutoffData(double cutoff, double relCutoff, double energy, std::vector<int> gridCounts)
    : _cutoff(cutoff), _relCutoff(relCutoff), _energy(energy), _gridCounts(std::move(gridCounts)){};

  // @brief Whether this class has the two given cutoffs
  bool has(double cutoff, double relCutoff) const {
    double eps = 1e-12;
    return std::fabs(cutoff - _cutoff) < eps && std::fabs(relCutoff - _relCutoff) < eps;
  };

  std::vector<int> getGridCounts() const {
    return _gridCounts;
  };

  double getEnergy() const {
    return _energy;
  };

 private:
  double _cutoff;
  double _relCutoff;
  double _energy;
  std::vector<int> _gridCounts;
};

/**
 * @class Cp2kGridDataContainer Cp2kGridData.h
 * @brief This class holds multiple Cp2kGridData objects
 */
class Cp2kCutoffDataContainer {
 public:
  Cp2kCutoffDataContainer() = default;

  void add(const Cp2kCutoffData& data) {
    _data.push_back(data);
  };

  // @brief Whether this container holds data with the two given cutoffs
  bool has(double cutoff, double relCutoff) const {
    return std::any_of(_data.begin(), _data.end(),
                       [&, cutoff, relCutoff](const auto& data) { return data.has(cutoff, relCutoff); });
  };

  Cp2kCutoffData getData(double cutoff, double relCutoff) const {
    for (const auto& data : _data) {
      if (data.has(cutoff, relCutoff)) {
        return data;
      }
    }
    throw std::runtime_error("Data for cutoffs " + std::to_string(cutoff) + " and " + std::to_string(relCutoff) +
                             " is not present in this container.");
  };

 private:
  std::vector<Cp2kCutoffData> _data;
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_EXTERNALQC_CP2KCUTOFFDATA_H
