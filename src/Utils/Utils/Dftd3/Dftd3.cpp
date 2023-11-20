/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Dftd3.h"
#include "Dftd3ReferencePairs.h"
#include <Utils/Math/AutomaticDifferentiation/MethodsHelpers.h>

namespace Scine {
namespace Utils {
namespace Dftd3 {

// Initialization of a D3 calculation.
void Dftd3::initialize(const AtomCollection& atomCollection, double s6, double s8, double dampingParam1,
                       double dampingParam2, Damping damping) {
  // initialize energy and gradients and the AtomicSecondDerivativesCollection, which holds
  // the gradients during the calculation.
  energy_ = 0.0;
  gradients_.resize(atomCollection.size(), 3);
  gradients_.setZero();
  UpToSecondDerivatives_ = AtomicSecondDerivativeCollection(atomCollection.size());
  UpToSecondDerivatives_.setZero();
  damping_ = damping;
  // Initialize D3 Parameters.
  parameters_ = Dftd3Parameters();
  parameters_.setS6(s6);
  parameters_.setS8(s8);
  if (damping_ == BJ) {
    parameters_.setA1(dampingParam1);
    parameters_.setA2(dampingParam2);
  }
  else if (damping_ == Zero) {
    parameters_.setSr(dampingParam1);
    parameters_.setA(dampingParam2);
  }
  else {
    throw std::runtime_error("Damping function not supported");
  }
  structure_.clear();

  // This loop converts the AtomCollection into a vector of Dftd3Atoms.
  for (int index = 0; index < atomCollection.size(); index++) {
    Dftd3Atom newAtom(atomCollection.getElement(index), atomCollection.getPosition(index));
    newAtom.setIndex(index);
    structure_.push_back(newAtom);
  }
}

// Getter for one cutoff radius R0 for one atom pair (atom 1 and atom 2)
double Dftd3::getR0(int atom1Index, int atom2Index) {
  return R0_(atom1Index, atom2Index);
}

// Getter for Matrix of C6 coefficients
Eigen::MatrixXd Dftd3::getC6Matrix() const {
  return C6_;
}

// Getter for the energy
double Dftd3::getEnergy() const {
  return energy_;
}

// Getter for the structure. Note, that this is a vector of Dftd3Atoms, not an AtomCollection.
std::vector<Dftd3Atom> Dftd3::getStructure() {
  return structure_;
}

// Setter for the structure. Note, that this is a vector of Dftd3Atoms, not an AtomCollection.
void Dftd3::setStructure(std::vector<Dftd3Atom> structure) {
  structure_ = std::move(structure);
}

// Getter for the Gradients.
GradientCollection Dftd3::getGradients() {
  return gradients_;
}

// This function calculates all the C6, C8 and R0 values and stores them in Eigen matrices.
void Dftd3::calculateValuesForConstants() {
  C6_.resize(structure_.size(), structure_.size());
  C8_.resize(structure_.size(), structure_.size());
  R0_.resize(structure_.size(), structure_.size());

  // Loop over all atom pairs only once
  for (const auto& atom1 : structure_) {
    for (const auto& atom2 : structure_) {
      if (atom1.getIndex() <= atom2.getIndex()) {
        continue;
      }

      // Calculate the values
      double c6 = calculateC6Coefficient(atom1, atom2);
      double c8 = calculateC8Coefficient(atom1, atom2, c6);
      double r0;

      if (damping_ == BJ) {
        r0 = std::sqrt(c8 / c6);
      }
      else if (damping_ == Zero) {
        r0 = parameters_.getR0Zero(atom1.getElementType(), atom2.getElementType());
      }
      else {
        throw std::runtime_error("Damping function not supported");
      }

      // Store the values into the Eigen matrices. Note, that C6_(i,j) = C6_(j,i).
      C6_(atom1.getIndex(), atom2.getIndex()) = c6;
      C6_(atom2.getIndex(), atom1.getIndex()) = c6;
      C8_(atom1.getIndex(), atom2.getIndex()) = c8;
      C8_(atom2.getIndex(), atom1.getIndex()) = c8;
      R0_(atom1.getIndex(), atom2.getIndex()) = r0;
      R0_(atom2.getIndex(), atom1.getIndex()) = r0;
    }
  }
}

// Calculation of coordination numbers for each atom. Those are fractional values in the Grimme D3 scheme.
double Dftd3::calculateCoordinationNumber(const Dftd3Atom& centralAtom) {
  double coordinationNumber = 0.0;
  // Get some parameters for the Dftd3 parameters.
  double k1 = parameters_.getK1();
  double covalentRadiusA = parameters_.getCovalentRadius(centralAtom.getElementType());

  // Loop over all atoms except the central atom and add their contribution to the
  // coordination number of the central atom.
  for (auto const& atom : structure_) {
    if (centralAtom.getIndex() == atom.getIndex()) {
      continue;
    }
    double distance = (atom.getPosition() - centralAtom.getPosition()).norm();
    double covalentRadiusB = parameters_.getCovalentRadius(atom.getElementType());

    // Formula for the contribution of one atom to the coord. numb. of the central atom.
    coordinationNumber += 1.0 / (1.0 + std::exp(-k1 * ((covalentRadiusA + covalentRadiusB) / distance - 1.0)));
  }

  return coordinationNumber;
}

// Derivative of the coord. number of atom 1 w.r.t. the distance between atom 1 and atom 2
double Dftd3::evaluateGradientOfCoordNumberWrtDistance(Dftd3Atom& atom1, Dftd3Atom& atom2) {
  // Make coordination number a Second1D object from the automatic differentiation library.
  AutomaticDifferentiation::Second1D coordinationNumber(0, 0, 0);
  double k1 = parameters_.getK1();
  double covalentRadiusA = parameters_.getCovalentRadius(atom1.getElementType());
  double covalentRadiusB = parameters_.getCovalentRadius(atom2.getElementType());

  // Make the distance between atom 1 and atom 2 a Second1D object.
  AutomaticDifferentiation::Second1D distance((atom2.getPosition() - atom1.getPosition()).norm(), 1, 0);

  // Almost all parts of the sum in the formula for the coord. numb. become zero. Only this term survives.
  coordinationNumber = 1.0 / (1.0 + exp(-k1 * ((covalentRadiusA + covalentRadiusB) / distance - 1.0)));

  // Return only the first derivative stored in the Second1D object.
  return coordinationNumber.first();
}

// Calculation of the C6 coefficient for the atom pair atom 1 and atom 2.
double Dftd3::calculateC6Coefficient(const Dftd3Atom& atom1, const Dftd3Atom& atom2) {
  double k3 = parameters_.getK3();

  // calculate W and Z
  double w = 0.0;
  double z = 0.0;

  // Get reference Pairs (3 x 75 dimensional std::array) for this element pair.
  const auto& referencePairs = parameters_.getReferencePairs(atom1.getElementType(), atom2.getElementType());

  // referencePairs[i][0] = Coordination Number 1
  // referencePairs[i][1] = Coordination Number 2
  // referencePairs[i][2] = C6 Coefficient
  for (int i = 0; i < 25; i++) {
    // Stop the for loop if one gets to the placeholder values.
    if (referencePairs[i][0] == Dftd3ReferencePairs::placeholder) {
      break;
    }
    // Implementation of the formula presented in the reference mentioned in the header file.
    double a = atom1.getCoordinationNumber() - referencePairs[i][0];
    double b = atom2.getCoordinationNumber() - referencePairs[i][1];
    double l = std::exp(-k3 * (a * a + b * b));
    w += l;
    z += referencePairs[i][2] * l;
  }

  double c6 = z / w;

  return c6;
}

// Derivative of the C6 coefficient for the atom pair atom 1 and atom 2 w.r.t. the coordination number of atom 1.
// The implementation is the same as in the "calculateC6Coefficient" function, but the values w, z and a
// are represented by a Second1D object. Thus, the C6 coefficient is obtained as a Second1D object.
// Only the first derivative part is returned.
double Dftd3::evaluateGradientOfC6WrtCoordNumber(Dftd3Atom& atom1, Dftd3Atom& atom2) {
  double k3 = parameters_.getK3();

  AutomaticDifferentiation::Second1D w(0, 0, 0);
  AutomaticDifferentiation::Second1D z(0, 0, 0);

  const auto& referencePairs = parameters_.getReferencePairs(atom1.getElementType(), atom2.getElementType());

  // referencePairs[i][0] = Coordination Number 1
  // referencePairs[i][1] = Coordination Number 2
  // referencePairs[i][2] = C6 Coefficient
  for (int i = 0; i < 25; i++) {
    if (referencePairs[i][0] == Dftd3ReferencePairs::placeholder) {
      break;
    }
    AutomaticDifferentiation::Second1D a(atom1.getCoordinationNumber() - referencePairs[i][0], 1, 0);
    double b = atom2.getCoordinationNumber() - referencePairs[i][1];
    auto l = exp(-k3 * (a * a + b * b));
    w += l;
    z += referencePairs[i][2] * l;
  }

  auto c6 = z / w;

  // Return first derivative.
  return c6.first();
}

// Calculation of the C8 coefficient for the atom pair atom 1 and atom 2. The C6 coefficient
// is needed as an argument.
double Dftd3::calculateC8Coefficient(const Dftd3Atom& atom1, const Dftd3Atom& atom2, double c6) {
  // Get some parameters for the element types of atom 1 and atom 2.
  double r2r41 = parameters_.getR2r4(atom1.getElementType());
  double r2r42 = parameters_.getR2r4(atom2.getElementType());

  // C8 is just a scaled version of C6.
  return 3.0 * c6 * r2r41 * r2r42;
}

// Evaluation of the D3 energy for just one atom pair (atom1, atom2).
double Dftd3::evaluateEnergy(Dftd3Atom& atom1, Dftd3Atom& atom2) {
  // Get the C6 and C8 coefficients as well as the cutoff radius, since they were pre-calculated for each atom pair.
  auto c6 = C6_(atom1.getIndex(), atom2.getIndex());
  auto c8 = C8_(atom1.getIndex(), atom2.getIndex());
  auto r0 = getR0(atom1.getIndex(), atom2.getIndex());
  // Get the scaling parameters
  auto s6 = parameters_.getS6();
  auto s8 = parameters_.getS8();

  // Calculate the distance between atom 1 and atom 2.
  double r = (atom1.getPosition() - atom2.getPosition()).norm();

  // Calculate the r^6 and r^8 terms of the energy expression separately.
  double damping6 = 0.0;
  double damping8 = 0.0;
  if (damping_ == BJ) {
    damping6 = evaluateBjDamping(r, r0, 6);
    damping8 = evaluateBjDamping(r, r0, 8);
  }
  else if (damping_ == Zero) {
    damping6 = evaluateZeroDamping(r, r0, 6);
    damping8 = evaluateZeroDamping(r, r0, 8);
  }
  else {
    throw std::runtime_error("Damping function not supported");
  }

  double e6 = s6 * damping6 * c6 / std::pow(r, 6);
  double e8 = s8 * damping8 * c8 / std::pow(r, 8);

  // Final step of the energy evaluation and return.
  return -(e6 + e8);
}

// This function evaluates the derivative of the energy w.r.t. the C6 coefficient.
// This is just dE/dC6 = E/C6.
double Dftd3::evaluateEnergyDividedByC6(Dftd3Atom& atom1, Dftd3Atom& atom2) {
  auto c6 = C6_(atom1.getIndex(), atom2.getIndex());
  // Use the "evaluateEnergy" function.
  return evaluateEnergy(atom1, atom2) / c6;
}

// Evaluate the partial derivative of the D3 energy for one atom pair w.r.t. their distance.
double Dftd3::evaluateGradientsWrtDistances(Dftd3Atom& atom1, Dftd3Atom& atom2) {
  // Get the C6 and C8 coefficients as well as the cutoff radius, since they were pre-calculated for each atom pair.
  auto c6 = C6_(atom1.getIndex(), atom2.getIndex());
  auto c8 = C8_(atom1.getIndex(), atom2.getIndex());
  auto r0 = getR0(atom1.getIndex(), atom2.getIndex());
  // Get the scaling parameters
  auto s6 = parameters_.getS6();
  auto s8 = parameters_.getS8();

  Displacement rVector = atom2.getPosition() - atom1.getPosition();

  // Distance is handled as a Second1D object.
  AutomaticDifferentiation::Second1D r(rVector.norm(), 1, 0);

  auto r2 = r * r;
  auto r3 = r * r2;
  auto r6 = r3 * r3;
  auto r8 = r6 * r2;

  // Calculate the r^6 and r^8 terms of the energy expression separately.
  auto damping6 = AutomaticDifferentiation::Second1D(0, 1, 0);
  auto damping8 = AutomaticDifferentiation::Second1D(0, 1, 0);
  if (damping_ == BJ) {
    damping6 = evaluateBjDamping(r, r0, 6);
    damping8 = evaluateBjDamping(r, r0, 8);
  }
  else if (damping_ == Zero) {
    damping6 = evaluateZeroDamping(r, r0, 6);
    damping8 = evaluateZeroDamping(r, r0, 8);
  }
  auto e6 = s6 * damping6 * c6 / r6;
  auto e8 = s8 * damping8 * c8 / r8;

  // Obtain the energy as a Second1D object.
  AutomaticDifferentiation::Second1D final = -(e6 + e8);

  // Return only the first derivative.
  return final.first();
}

// Evaluate the contribution of one atom pair to the total gradients.
void Dftd3::evaluateGradients(Dftd3Atom& atom1, Dftd3Atom& atom2, double& dE_drij, std::vector<double>& dc6,
                              AtomicSecondDerivativeCollection& UpToSecondDerivatives) {
  // Set up a Second1D object for the energy. The zeroth and second derivative are not set, because one is only
  // interested in the first derivative. The expression for the first derivative is set up by hand using the
  // terms dE/drij and dc6 (meanings of those values are described in the "calculate" function).
  AutomaticDifferentiation::Second1D energySecond1D(
      0, dE_drij + dCN_drij_(atom1.getIndex(), atom2.getIndex()) * (dc6[atom1.getIndex()] + dc6[atom2.getIndex()]), 0);

  // Define difference vector.
  auto rVector = atom2.getPosition() - atom1.getPosition();
  // Get a Second3D object from the Second1D object under consideration of the difference vector.
  auto v = AutomaticDifferentiation::get3Dfrom1D<DerivativeOrder::Two>(energySecond1D, rVector);

  // Transfer the gradient contributions for atom 1 and atom 2 into the UpToSecondDerivatives container.
  UpToSecondDerivatives[atom2.getIndex()] += v;
  UpToSecondDerivatives[atom1.getIndex()] += v.opposite();
}

// This function does a D3 calculation. The derivative type (none or first) is given as an argument.
void Dftd3::calculate(Derivative d) {
  // Calculate and set coordination numbers for all atoms
  for (auto& atom : structure_) {
    double coordinationNumber = calculateCoordinationNumber(atom);
    atom.setCoordinationNumber(coordinationNumber);
  }

  // Calculate values for C6, C8 and R0.
  calculateValuesForConstants();

  // Set matrices dC6/dCN and dCN/drij.
  // (dC6/dCN)ij = dC6(ij)/dCN(i)
  // (dCN/drij)ij = dCN(i)/drij
  if (d != Derivative::None) {
    dC6_dCN_.resize(structure_.size(), structure_.size());
    dCN_drij_.resize(structure_.size(), structure_.size());
    // Loop over all atom pairs.
    for (auto& atom1 : structure_) {
      for (auto& atom2 : structure_) {
        if (atom1.getIndex() == atom2.getIndex()) {
          continue;
        }

        dC6_dCN_(atom1.getIndex(), atom2.getIndex()) = evaluateGradientOfC6WrtCoordNumber(atom1, atom2);
        dCN_drij_(atom1.getIndex(), atom2.getIndex()) = evaluateGradientOfCoordNumberWrtDistance(atom1, atom2);
      }
    }
  }

  // Set dc6 values for each atom.
  // Definition of a dc6 value:
  // dc6(i) = sum over all (j != i) of dC6(ij)/dCN(i) * dE(ij)/dC6(ij)
  std::vector<double> dc6(structure_.size(), 0.0);
  if (d != Derivative::None) {
    // Loop over all atom pairs only once
    for (auto& atom1 : structure_) {
      for (auto& atom2 : structure_) {
        if (atom1.getIndex() >= atom2.getIndex()) {
          continue;
        }

        dc6[atom1.getIndex()] += dC6_dCN_(atom1.getIndex(), atom2.getIndex()) * evaluateEnergyDividedByC6(atom1, atom2);
        dc6[atom2.getIndex()] += dC6_dCN_(atom2.getIndex(), atom1.getIndex()) * evaluateEnergyDividedByC6(atom1, atom2);
      }
    }
  }

  // Loop over all atom pairs once to calculate the energy and optionally derivatives.
  for (auto& atom1 : structure_) {
    for (auto& atom2 : structure_) {
      if (atom1.getIndex() >= atom2.getIndex()) {
        continue;
      }
      // Calculate the energy for this atom pair.
      energy_ += evaluateEnergy(atom1, atom2);
      // If a first derivative is required, calculate it.
      if (d != Derivative::None) {
        // Calculate the partial derivatives of the energy w.r.t. the distance for one atom pair.
        double dE_drij = evaluateGradientsWrtDistances(atom1, atom2);
        // Evaluate the contribution of the atom pair atom 1 and atom 2 to the total gradients.
        // The UpToSecondDerivatives object is given by reference and modified within the function.
        evaluateGradients(atom1, atom2, dE_drij, dc6, UpToSecondDerivatives_);
      }
    }
  }

  // Convert the UpToSecondDerivatives matrix to the GradientCollection once it is complete.
  if (d != Derivative::None) {
    const unsigned N = structure_.size();
#pragma omp parallel for
    for (int i = 0; i < N; ++i) {
      gradients_.row(i) = Gradient(UpToSecondDerivatives_[i].deriv());
    }
  }
}

} // namespace Dftd3
} // namespace Utils
} // namespace Scine
