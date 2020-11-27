=========
Changelog
=========

Release 3.0.1
=============

- Update development-utils

Release 3.0.0
=============

- Enable ``conan`` builds and PyPI releases
- Add Newton-Trajectory type optimizer
- Add statistics and machine learning tools (PCA, k-fold cross-validation, kernel ridge regression)
- Add chemical representations for machine-learned force fields
- Charge Model 5 corrections for Hirshfeld atomic partial charges
- Add various conceptual DFT quantities
- Enable access to the density matrix and GTOs in Python
- Improve Dimer transition state search algorithm
- Add implicit solvation option to ORCA and Gaussian interfaces
- Add python bindings sphinx documentation
- Separate ``Settings`` from its base ``ValueCollection`` in python bindings

Release 2.0.0
=============

- Add support for internal coordinates
- Add interface to Gaussian
- Improve ORCA interface (and make compatible with ORCA 4.2.0)
- Add BFGS optimizer and G-DIIS convergence accelerator
- Improve Bofill transition state search algorithm
- Various bugfixes and improvements

Release 1.0.1
=============

Hotfix to allow compilation on OSX using Clang.

Release 1.0.0
=============

Initial release with all necessary functionality to support Sparrow and ReaDuct.
Among other things, this includes:

- Analytic evaluation of gradients
- Calculation of bond orders
- Interface to the ORCA quantum chemistry program
- Numerical Hessian calculator
- Optimizers to find minima and transition states on the PES
- Python bindings
- SCF algorithm (including convergence accelerators such as DIIS)
