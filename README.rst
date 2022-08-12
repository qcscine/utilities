SCINE - Utilities
=================

Introduction
------------

This repository contains functionality which is used in most SCINE modules.
It is vital for the correct functioning of all SCINE modules, but it does not
directly provide any functionality for end users. Therefore, only developers
of SCINE need to directly interact with this repository.

License and Copyright Information
---------------------------------

For license and copyright information, see the file ``LICENSE.txt`` in this
directory.

Installation and Usage
----------------------

As a user of one of the SCINE modules (such as Sparrow), you do not need
to set up SCINE Utilities yourself. Instead, this is done as part of the
installation process of the respective SCINE module. Therefore, the following
instructions are only necessary for developers of SCINE, or those interfacing
this library directly.

Dependencies
............

Required software, minimum required versions in brackets, for this SCINE project are:

- A C++ compiler supporting the C++17 standard
- CMake (3.9)
- Boost (1.65.0)
- Eigen3 (3.3.2)

Installation
............

In order to compile this as a SCINE developer, execute the following
commands::

    git submodule init
    git submodule update
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../inst ..
    make -j 4
    make UtilsOSDocumentation
    make test
    make install

Known Issues
------------

Due to a `bug <https://github.com/boostorg/process/issues/62>`_ in boost 
version 1.69.0, linkage against this version will result in a leakage
of file handles if the ORCA and Gaussian Calculators are called repeatedly.

Support and Contact
-------------------

In case you should encounter problems or bugs, please write a short message
to scine@phys.chem.ethz.ch.

Third-Party Libraries Used
--------------------------

SCINE Utilities makes use of the following third-party libraries:

- `Boost <https://www.boost.org/>`_
- `Eigen <http://eigen.tuxfamily.org>`_
- `Google Test <https://github.com/google/googletest>`_
- `IRC <https://github.com/rmeli/irc>`_
- `LBFGSB <https://github.com/yixuan/LBFGSpp>`_
- `pybind11 <https://github.com/pybind/pybind11>`_
- `yaml-cpp <https://github.com/jbeder/yaml-cpp>`_
