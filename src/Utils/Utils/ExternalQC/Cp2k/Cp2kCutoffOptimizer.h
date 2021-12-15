/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_CP2KCUTOFFOPTIMIZER_H
#define UTILS_CP2KCUTOFFOPTIMIZER_H

#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/ExternalQC/Cp2k/Cp2kCalculator.h>
#include <Utils/ExternalQC/Cp2k/Cp2kCutoffData.h>
#include <Utils/Settings.h>
#include <algorithm>
#include <numeric>

namespace Scine {
namespace Utils {
namespace ExternalQC {

class Cp2kCutoffOptimizer {
 public:
  Cp2kCutoffOptimizer(Core::Calculator& calculator, Cp2kCutoffDataContainer dataContainer = Cp2kCutoffDataContainer())
    : _calculator(calculator), _dataContainer(std::move(dataContainer)){};
  /*
   * @brief Determine the optimal planeWaveCutoff and relative multigrid cutoff for the given system
   *
   * Performs multiple single point calculations with different grid parameters with just a single SCF
   * cycle. It automatically determines the optimal grid parameters for this structure within the given margins and
   * writes them into the settings.
   *
   * @param energyEps             The threshold for sufficient energy convergence based on the grid cutoffs in hartree.
   * @param distributionEpsFactor The threshold for a sufficient distribution of the Gaussian functions on the
   *                              different grids given as a factor of the ideal distribution. The ideal distribution
   *                              would be that each sub grid has the exact same percentage of Gaussian of the total
   *                              number of Gaussian functions (1 / nGrids). This threshold determines that no subgrid
   *                              has a lower percentage of Gaussian functions than the threshold multiplied with the
   *                              ideal percentage (distributionEpsFactor / nGrids).
   * @param startCutoff           The cutoff with which the optimization will be started. It serves as a first reference.
   *                              If the cutoff right below that already deviates from this cutoff, the cutoff is
   *                              increased and therefore the reference is changed
   * @param startRelCutoff        Just like startCutoff, but for the relCutoff
   */
  void determineOptimalGridCutoffs(double energyEps = 1e-8, double distributionEpsFactor = 0.01,
                                   double startCutoff = 500.0, double startRelCutoff = 100.0);

 private:
  // @brief converge either the cutoff or relative cutoff to the energy threshold
  double convergeCutoff(double startCutoff, double fixedCutoff = 100.0, bool convergeRelCutoff = false);
  // @brief select cutoffs that also fulfill a proper distribution of Gaussian functions on the different grids
  std::pair<double, double> convergeDistribution(double cutoff, double relCutoff);
  // @brief retrieve data from container if available or else calculate it
  Cp2kCutoffData getGridData(double cutoff, double relCutoff);
  // @brief check against hard limit and throw exception
  void avoidInfiniteLoop(double maxCutoff, double hardLimit, double fixedCutoff, bool convergeRelCutoff) const;

  Core::Calculator& _calculator;
  Cp2kCutoffDataContainer _dataContainer;
  double _cutoffStepSize = 50.0;
  double _relCutoffStepSize = 10.0;
  double _cutoffHardLimit = 2000.0;
  double _relCutoffHardLimit = 500.0;
  double _energyEps;
  double _distributionEpsFactor;
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_CP2KCUTOFFOPTIMIZER_H
