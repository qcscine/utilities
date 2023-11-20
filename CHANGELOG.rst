Changelog
=========

All changes are grouped into one of three categories: those changes that are considered
being related to features, meaning that they add or update new functionalities,
those that are considered important technical updates, and all other notable updates.
It is intended that only the first two groups (``New Features and Feature Updates`` and
``Important Technical Changes``) are important for the average user, while
the last one is mainly aimed at developers and users that link deeply into the code.

Release 9.0.0
-------------

New Features and Feature Updates
................................
- Enable calculation of 57-Fe Mössbauer parameters with Orca.
- Enable DLPNO-CCSD(T)-F12 calculations with Orca.
- Added a MRCC calculator (requires MRCC executables at the `MRCC_BINARY_PATH` path
  environment variable).
  -  Supported methods: DFT, HF, MP2 [LNO-MP2], CC [LNO-CCSD, LNO-CCSD(T)]
- Added transition state optimizer for QM/MM
- Molecular trajectories in the PDB file format can now be read.

Important Technical Changes
...........................
- Improve support for compilation on Windows (MSVC)

Other
.....
- Update address in license

Release 8.0.0
-------------

New Features and Feature Updates
................................
- Turbomole calculator can now calculate Hessians numerically.
- The DFT grid for the Turbomole calculator can now be varied.
- The cavity construction for implicit solvation with the Turbomole calculator can now be modified by setting the points and segments per atom.
- Update Calculator interface helpers to accommodate pure Python3 Calculators.
- Added distortion of structure along a set of vibrational mode(s) and harmonic inversion point calculation.

Important Technical Changes
...........................
- Changed Python state definition of ValueCollection to avoid faulty state for empty ValueCollections.

Release 7.0.0
-------------

New Features and Feature Updates
................................
- A custom implicit solvent may be defined for the Turbomole input creator through
  its dielectric constant and probe radius.
- Turbomole calculator can now carry out Hartree-Fock calculations.
- Orca calculator can now carry out broken-symmetry DFT calculations.
- Improve comparisons of periodic data structures.
- Ease communication between forked processes for monitoring calculations.
- Various bugfixes and improvements
- Added data structures for interfacing with integral libraries.
- The BondDetectors can alternatively rely on Van der Waals (VdW) radii instead of covalent radii.
- The PeriodicSystem generates its bond orders now per default based on VdW radii between solid state atoms.
- Turbomole calculator can now calculate excited states.
- The PdbStreamHandler can read multiple models in a pdb file (encoded by the MODEL tag in the file).
- Add Python bindings for D3/D3BJ dispersion corrections.
- Add additional setting to external QC calculators that allows to enforce the set SCF convergence criterion.
- The pressure for free energy calculations may now be changed. Current default 1 atm.

Important Technical Changes
...........................
- Add Spglib dependency.
- Increase pickle support of Python bindings.
- Optimization related setting names are gathered in a separate namespace and removed from the optimizers.

Release 6.0.0
-------------

New Features and Feature Updates
................................
- Newton Trajectory: Added new extraction options and improved eta bonds in NT2
- Added more solvents for the Turbomole input creator
- Added the Pauling electronegativity scale. This is available through the ElementInfo.
- Added solvate function to place any mixture of solvents around solute.

Important Technical Changes
...........................
- PeriodicSystem canonicalizes the given PeriodicBoundaries to ensure read/write stability
- Calculators:
  - For the Turbomole calculator, allow the SCF damping value to be specified exactly (instead of predefined
    settings "default", "low", "medium", and "high")
  - Add "PBEH-3C" and "B97-3C" as supported method families for the ORCA calculator
- Settings for placing solvent molecules are summarized in ``SolventPlacementSettings``, all solvate functions
  take only this structure as input. This change is not backwards compatible.

Development Changes
...................
- Code deduplication

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
