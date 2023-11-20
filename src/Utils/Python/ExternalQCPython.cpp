/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/ExternalQC/Cp2k/Cp2kCutoffOptimizer.h>
#include <pybind11/pybind11.h>

using namespace Scine::Utils::ExternalQC;

void init_cp2k_external_qc(pybind11::module& m) {
  pybind11::class_<Cp2kCutoffOptimizer> cutoff_optimizer(m, "Cp2kCutoffOptimizer");
  cutoff_optimizer.def(pybind11::init<Scine::Core::Calculator&>(), pybind11::arg("cp2k_calculator"),
                       "Initialize with a particular calculator with settings and a structure.");

  cutoff_optimizer.def("determine_optimal_grid_cutoffs", &Cp2kCutoffOptimizer::determineOptimalGridCutoffs,
                       pybind11::arg("energy_accuracy") = 1e-8, pybind11::arg("distribution_factor_accuracy") = 0.01,
                       pybind11::arg("start_cutoff") = 500.0, pybind11::arg("start_rel_cutoff") = 100.0,
                       R"delim(
                       Determine cutoff and relCutoff and set it directly in the settings of the calculator.
                       'energy_accuracy' is the threshold for sufficient energy convergence based on the grid cutoffs
                       in hartree.
                       'distribution_factor_accuracy' is the threshold for a sufficient distribution of the Gaussian
                       functions on the different grids given as a factor of the ideal distribution. The ideal
                       distribution would be that each sub grid has the exact same percentage of Gaussian of the total
                       number of Gaussian functions (1 / n_grids). This threshold determines that no subgrid
                       has a lower percentage of Gaussian functions than the threshold multiplied with the
                       ideal percentage (distributionEpsFactor / nGrids).
                       'start_cutoff' ist the cutoff with which the optimization will be started. It serves as a first
                       reference. If the cutoff right below that already deviates from this cutoff, the cutoff is
                       increased and therefore the reference is changed.
                       'start_rel_cutoff' is just like 'start_cutoff', but for the relCutoff.
                       )delim");
}
