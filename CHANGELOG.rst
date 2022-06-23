Changelog
=========

All changes are grouped into one of three categories: those changes that are considered
being related to features, meaning that they add or update new functionalities,
those that are considered important technical updates, and all other notable updates.
It is intended that only the first two groups (``New Features and Feature Updates`` and
``Important Technical Changes``) are important for the average user, while
the last one is mainly aimed at developers and users that link deeply into the code.

Release 5.0.0
-------------

New Features and Feature Updates
................................
- Add second version of Newton Trajectory algorithm (NT2) for elementary step induction
- Add Gaussian process regression
- Improve BFGS: guarantees positive-definite approximate Hessian
- B-Splines: interpolation and compression of trajectories including energy
- Enable optimization of a geometry together with its unit cell
  - CP2K:
    - Add xTB support
    - Add stress tensor
- Enable point charge embedding for Turbomole
- Add linear sum assignment algorithm
- Add functionality to evaluate the spin contamination of a single determinant
- Add Python bindings for WavefunctionOutputGenerator and casting utilities
  from Calculator and CalculatorWithReference.

Important Technical Changes
...........................
- Calculators:
  - Harmonize dispersion correction input

Development Changes
...................
- Code deduplication

Release 4.0.0
-------------

New Features and Feature Updates
................................
- Improve GDIIS: numerical stability
- Improve BFGS: add automatic damping
- Improve EVF and Bofill: allows to select mode and follows mode independent of order
- Add periodic boundary conditions
- Add support for CP2K
- Improve Gaussian interface: allows to reuse SCF results as guesses in subsequent calculations and
  to retrieve molecular orbital coefficients
- Add support for Turbomole
- Improve MD: Fix a bug regarding the time step size, check gradient
  calculations for SCF convergence, add the option to use bias potentials and
  add a stochastic dynamics integrator

Important Technical Changes
...........................
- Add Python bindings for CalculatorWithReference
- Add Log accessor to Python bindings of Calculator and CalculatorWithReference
- Add Python bindings for Davidson diagonalizer with possibility of having
  custom sigma vector evaluators/preconditioners
- Add functions to get all closest atoms within a certain distance and to
  build an atom pair-list
- Distinguish now between true internal, true Cartesian, and Cartesian with removed
  translation and rotation coordinate systems
- Add Python bindings for ThermochemistryCalculator to calculate thermodynamic properties from a Hessian in Python
- Add Python bindings for SettingsNames
- Add support for the SMD solvation model in ORCA
- Add the option to obtain gradients from a CalculatorWithReference in MD simulations

Development Changes
...................
- Refactoring of GeometryUtilities into sub-namespaces
- Add data structures needed for downstream methods that are general to linear response methods
- Remove Logger option for downstream LcaoMethods as it can be accessed through the calculator interface
- Refactor Davidson diagonalizers:

  - Create IterativeDiagonalizer interface
  - Create KrylovDiagonalizer interface, inheriting from IterativeDiagonalizer
  - Create the 2 versions, NonOrthogonalDavidson and OrthogonalDavidson
  - Add Python bindings for OrthogonalDavidson and NonOrthogonalDavidson,
    tested in Python and added an example on how to extend the SigmaVectorEvaluator
    to customize the Davidson directly in Python

Release 3.0.1
-------------

Important Technical Changes
...........................

- Update development-utils

Release 3.0.0
-------------

New Features and Feature Updates
................................
- Add Newton-Trajectory reaction search optimizer
- Improve Dimer transition state search algorithm
- Improve BFGS/GDIIS geometry optimization algorithm
- Add statistics and machine learning tools (PCA, k-fold cross-validation, kernel ridge regression)
- Add chemical representations for machine-learned force fields
- Add possibility to generate Charge Model 5 (CM5) corrections for Hirshfeld atomic partial charges
- Add various conceptual DFT quantities
- Add implicit solvation options to ORCA and Gaussian interfaces

Important Technical Changes
...........................
- Enable ``conan`` builds and PyPI releases
- Add Python bindings sphinx documentation

Development Changes
...................
- Enable access to the density matrix and GTOs in Python
- Separate ``Settings`` from its base ``ValueCollection`` in Python bindings
- Add Python bindings for molecular dynamics simulations
- Rework the Python wrapper for ``Settings``, ``ValueCollection`` and ``DescriptorCollection``
- Add a ``TestCalculator`` and module that implements the ``Test`` method to allow mocked calls
  to QC programs. (Uses a modified Lennard-Jones potential)

Release 2.0.0
-------------

- Add support for internal coordinates
- Add interface to Gaussian
- Improve ORCA interface (and make compatible with ORCA 4.2.0)
- Add BFGS optimizer and G-DIIS convergence accelerator
- Improve Bofill transition state search algorithm
- Various bugfixes and improvements

Release 1.0.1
-------------

Hotfix to allow compilation on OSX using Clang.

Release 1.0.0
-------------

Initial release with all necessary functionality to support Sparrow and ReaDuct.
Among other things, this includes:

- Analytic evaluation of gradients
- Calculation of bond orders
- Interface to the ORCA quantum chemistry program
- Numerical Hessian calculator
- Optimizers to find minima and transition states on the PES
- Python bindings
- SCF algorithm (including convergence accelerators such as DIIS)
