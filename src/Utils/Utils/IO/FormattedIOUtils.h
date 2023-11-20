/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_FORMATTEDIOUTILS_H
#define UTILS_FORMATTEDIOUTILS_H

#include <Utils/DataStructures/AtomsOrbitalsIndexes.h>
#include <Utils/GeometricDerivatives/NormalModesContainer.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/Geometry/ElementTypes.h>
#include <Utils/Math/IterativeDiagonalizer/SpinAdaptedEigenContainer.h>
#include <Utils/Properties/Thermochemistry/ThermochemistryCalculator.h>
#include <Utils/TimeDependent/TransitionDipoleCalculator.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <fstream>
#include <iomanip>
#include <ostream>

namespace Scine {
namespace Utils {

template<class MatrixType>
inline std::pair<int, int> getMatrixDimensions(const MatrixType& matrix) {
  return {matrix.size(), matrix[0].size()};
}

template<>
inline std::pair<int, int> getMatrixDimensions(const Eigen::MatrixXd& matrix) {
  return {matrix.rows(), matrix.cols()};
}

template<>
inline std::pair<int, int> getMatrixDimensions(const Eigen::Matrix3d& /* matrix */) {
  return {3, 3};
}

template<>
inline std::pair<int, int> getMatrixDimensions(const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& matrix) {
  return {matrix.rows(), matrix.cols()};
}

template<class MatrixType>
inline double getElement(int row, int col, const MatrixType& matrix) {
  return matrix[row][col];
}

template<>
inline double getElement(int row, int col, const Eigen::Matrix3d& matrix) {
  return matrix.col(col)(row);
}

template<>
inline double getElement(int row, int col, const Eigen::MatrixXd& matrix) {
  return matrix.col(col)(row);
}

template<>
inline double getElement(int row, int col, const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& matrix) {
  return matrix.col(col)(row);
}

template<>
inline double getElement(int row, int col, const NormalModesContainer& matrix) {
  auto displCollection = matrix.getMode(col);
  return displCollection.row(row / 3)(row % 3);
}

template<class MatrixType>
inline void matrixPrettyPrint(std::ostream& out, const MatrixType& matrix) {
  auto dimension = getMatrixDimensions(matrix);
  auto rest = dimension.second % 5;
  auto columnsToRest = dimension.second - rest;

  out << std::scientific << std::right;
  // Print the matrix in blocks of 5 columns
  for (int col = 0; col < dimension.second - 5; col += 5) {
    // print header line
    out << std::setw(10) << " ";
    for (int i = col; i < col + 5; ++i) {
      out << std::setw(15) << i + 1;
    }
    out << std::endl;
    // Print the whole columns
    for (int row = 0; row < dimension.first; ++row) {
      out << std::setw(10) << row + 1;
      for (int i = col; i < col + 5; ++i) {
        out << std::setw(15) << getElement(row, i, matrix);
      }
      out << std::endl;
    }
    out << std::endl;
  }

  // print rest of columns
  if (rest != 0) {
    // print header line
    out << std::setw(10) << " ";
    for (auto i = 0; i < rest; ++i) {
      out << std::setw(15) << columnsToRest + i + 1;
    }
    out << std::endl;
    for (int row = 0; row < dimension.first; ++row) {
      out << std::setw(10) << row + 1;
      for (int i = 0; i < rest; ++i) {
        out << std::setw(15) << getElement(row, columnsToRest + i, matrix);
      }
      out << std::endl;
    }
    out << std::defaultfloat << std::endl;
  }
}

inline void printElement(std::ostream& out, const Eigen::SparseMatrix<double>::InnerIterator& it) {
  out << std::setw(7) << std::left << "<" + std::to_string(it.row()) + ">";
  out << std::setw(7) << std::left << "<" + std::to_string(it.col()) + ">";
  out << std::setw(7) << std::left << std::to_string(it.value());
  out << std::setw(4) << std::left << " ";
  out << std::right;
}

inline void matrixPrettyPrint(std::ostream& out, const Eigen::SparseMatrix<double>& matrix, double threshold) {
  int nEl = 0;
  int nRows = 0;
  int headerCount = 0;
  const int nCols = 5;
  // print header line
  out << std::setw(10) << std::right << std::scientific << " ";
  for (++headerCount; headerCount <= nCols; ++headerCount) {
    out << std::setw(26) << std::left << headerCount;
  }
  out << std::endl;
  for (int k = 0; k < matrix.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(matrix, k); it; ++it) {
      if (it.value() > threshold) {
        if ((nEl % nCols) == 0) {
          out << std::endl;
          out << std::setw(10) << ++nRows;
        }
        if (nRows == 50) { // print line header
          nRows = 0;
          out << std::endl;
          out << std::endl;
          // print header line
          out << std::setw(10) << " ";
          for (++headerCount; (headerCount % 5) != 0; ++headerCount) {
            out << std::setw(26) << std::left << headerCount;
          }
          out << std::endl;
        }
        printElement(out, it);
        ++nEl;
      }
    }
  }
  out << std::defaultfloat;
}

inline void molecularOrbitalsPrettyPrint(std::ostream& out, const Eigen::MatrixXd& molecularOrbitals,
                                         const std::vector<double>& orbitalEnergies, const AtomsOrbitalsIndexes& aoInfo,
                                         const ElementTypeCollection& elementTypes) {
  auto dimension = molecularOrbitals.cols();
  auto rest = dimension % 5;
  auto columnsToRest = dimension - rest;
  // s, px, py, pz, dz2, dxz, dyz, dx2-y2, dxy
  std::vector<std::string> ao_symbols = {"s ", "px", "py", "pz", "dz2", "dxz", "dyz", "dx2-y2", "dxy"};
  out << std::endl
      << "======================================================================================" << std::endl
      << "Molecular Orbitals:" << std::endl
      << std::endl;
  out << std::scientific << std::right;
  // Print the matrix in blocks of 5 columns
  int col = 0;
  int colsToPrint = 5;
  while (col < dimension) {
    if (col >= columnsToRest) {
      colsToPrint = rest;
    }
    // print header line
    out << std::setw(10) << " ";
    for (int i = col; i < col + colsToPrint; ++i) {
      out << std::setw(15) << i;
    }
    // Print energies
    out << std::endl;
    out << std::setw(10) << "Energy:";
    for (int i = col; i < col + colsToPrint; ++i) {
      out << std::setw(15) << std::defaultfloat << orbitalEnergies[i];
    }
    out << std::endl << std::endl << std::fixed;
    // Print the whole columns
    int nAtoms = aoInfo.getNAtoms();
    int row = 0;
    for (int atom = 0; atom < nAtoms; atom++) {
      int nAOs = aoInfo.getNOrbitals(atom);
      std::string rowName = std::to_string(atom) + ElementInfo::symbol(elementTypes[atom]);
      for (int ao = 0; ao < nAOs; ao++) {
        out << std::setw(10) << rowName + " " + ao_symbols[ao];
        for (int i = col; i < col + colsToPrint; ++i) {
          out << std::setw(15) << molecularOrbitals(row, i);
        }
        row++;
        out << std::endl;
      }
    }
    out << std::endl;
    col += colsToPrint;
  }
  out << "======================================================================================" << std::endl
      << std::endl;
  out << std::defaultfloat;
}

inline void matrixPrettyPrint(std::ostream& out, const NormalModesContainer& container,
                              const ElementTypeCollection& elementTypes) {
  auto dimension = container.size();
  auto rest = dimension % 5;
  auto rows = container.getMode(0).size();
  auto columnsToRest = dimension - rest;
  auto waveNumbers = container.getWaveNumbers();

  out << std::scientific << std::right;
  // Print the matrix in blocks of 5 columns
  for (int col = 0; col < columnsToRest; col += 5) {
    // print header line
    out << std::setw(10) << " ";
    for (int i = col; i < col + 5; ++i) {
      out << std::setw(15) << i + 1;
    }
    out << std::endl;
    out << std::endl;
    // Print wavenumbers
    out << std::setw(10) << "Frequency:";
    for (int i = col; i < col + 5; ++i) {
      out << std::setw(15) << std::defaultfloat << waveNumbers[i];
    }
    out << std::endl << std::endl << std::scientific;
    // Print the whole columns
    for (int row = 0; row < rows; ++row) {
      std::string rowName = std::to_string(row / 3 + 1) + " " + ElementInfo::symbol(elementTypes[row / 3]);

      if (row % 3 == 0) {
        rowName += " x";
      }
      else if (row % 3 == 1) {
        rowName += " y";
      }
      else {
        rowName += " z";
      }
      out << std::setw(10) << rowName;
      for (int i = col; i < col + 5; ++i) {
        out << std::setw(15) << getElement(row, i, container);
      }
      out << std::endl;
    }
    out << std::endl;
  }

  // print rest of columns
  if (rest != 0) {
    // print header line
    out << std::setw(10) << " ";
    for (auto i = 0; i < rest; ++i) {
      out << std::setw(15) << columnsToRest + i + 1;
    }
    out << std::endl << std::endl;

    // Print wavenumbers
    out << std::setw(10) << "Frequency:";
    for (int i = 0; i < rest; ++i) {
      out << std::setw(15) << std::defaultfloat << waveNumbers[columnsToRest + i];
    }
    out << std::endl << std::endl;
    for (int row = 0; row < rows; ++row) {
      std::string rowName = std::to_string(row / 3 + 1) + " " + ElementInfo::symbol(elementTypes[row / 3]);
      if (row % 3 == 0) {
        rowName += " x";
      }
      else if (row % 3 == 1) {
        rowName += " y";
      }
      else {
        rowName += " z";
      }
      out << std::setw(10) << rowName;
      for (int i = 0; i < rest; ++i) {
        out << std::setw(15) << std::scientific << getElement(row, columnsToRest + i, container);
      }
      out << std::endl;
    }
    out << std::defaultfloat << std::endl;
  }
}

inline void prettyPrint(std::ostream& out, const ThermochemicalComponentsContainer& thermochemicalComponents) {
  // header line
  out << std::setw(20) << " " << std::setw(27) << std::left << "Molecular symmetry number:" << std::setw(10)
      << thermochemicalComponents.overall.symmetryNumber << std::endl;
  out << std::scientific << std::endl;
  out << std::setw(20) << "Component:" << std::setw(35) << "ZPVE [kJ/mol]" << std::setw(35) << "Enthalpy [kJ/mol]"
      << std::setw(35) << "Entropy [kJ/(K mol)]" << std::endl;
  // Components
  out << std::setw(20) << "Vibration:" << std::setw(35)
      << thermochemicalComponents.vibrationalComponent.zeroPointVibrationalEnergy << std::setw(35)
      << thermochemicalComponents.vibrationalComponent.enthalpy << std::setw(35)
      << thermochemicalComponents.vibrationalComponent.entropy << std::endl;
  out << std::setw(20) << "Rotation:" << std::setw(35) << thermochemicalComponents.rotationalComponent.zeroPointVibrationalEnergy
      << std::setw(35) << thermochemicalComponents.rotationalComponent.enthalpy << std::setw(35)
      << thermochemicalComponents.rotationalComponent.entropy << std::endl;
  out << std::setw(20) << "Translation:" << std::setw(35)
      << thermochemicalComponents.translationalComponent.zeroPointVibrationalEnergy << std::setw(35)
      << thermochemicalComponents.translationalComponent.enthalpy << std::setw(35)
      << thermochemicalComponents.translationalComponent.entropy << std::endl;
  out << std::setw(20) << "Electronic:" << std::setw(35)
      << thermochemicalComponents.electronicComponent.zeroPointVibrationalEnergy << std::setw(35)
      << thermochemicalComponents.electronicComponent.enthalpy << std::setw(35)
      << thermochemicalComponents.electronicComponent.entropy << std::endl;
  out << std::setw(20) << "Overall:" << std::setw(35) << thermochemicalComponents.overall.zeroPointVibrationalEnergy
      << std::setw(35) << thermochemicalComponents.overall.enthalpy << std::setw(35)
      << thermochemicalComponents.overall.entropy << std::endl;
  out << std::endl;
  out << std::setw(20) << " " << std::setw(35) << "Heat Capacity Cp [kJ/(K mol)]" << std::setw(35)
      << "Heat Capacity Cv [kJ/(K mol)]" << std::setw(35) << "Gibbs Free Energy [kJ/mol]" << std::endl;
  // Components
  out << std::setw(20) << "Vibration:" << std::setw(35) << thermochemicalComponents.vibrationalComponent.heatCapacityP
      << std::setw(35) << thermochemicalComponents.vibrationalComponent.heatCapacityV << std::setw(35)
      << thermochemicalComponents.vibrationalComponent.gibbsFreeEnergy << std::endl;
  out << std::setw(20) << "Rotation:" << std::setw(35) << thermochemicalComponents.rotationalComponent.heatCapacityP
      << std::setw(35) << thermochemicalComponents.rotationalComponent.heatCapacityV << std::setw(35)
      << thermochemicalComponents.rotationalComponent.gibbsFreeEnergy << std::endl;
  out << std::setw(20) << "Translation:" << std::setw(35) << thermochemicalComponents.translationalComponent.heatCapacityP
      << std::setw(35) << thermochemicalComponents.translationalComponent.heatCapacityV << std::setw(35)
      << thermochemicalComponents.translationalComponent.gibbsFreeEnergy << std::endl;
  out << std::setw(20) << "Electronic:" << std::setw(35) << thermochemicalComponents.electronicComponent.heatCapacityP
      << std::setw(35) << thermochemicalComponents.electronicComponent.heatCapacityV << std::setw(35)
      << thermochemicalComponents.electronicComponent.gibbsFreeEnergy << std::endl;
  out << std::setw(20) << "Overall:" << std::setw(35) << thermochemicalComponents.overall.heatCapacityP << std::setw(35)
      << thermochemicalComponents.overall.heatCapacityV << std::setw(35)
      << thermochemicalComponents.overall.gibbsFreeEnergy << std::endl;
  out << std::endl;
}

inline void prettyPrint(std::ostream& out, const SpinAdaptedElectronicTransitionResult& result) {
  // Lambda to avoid duplication for singlet and triplet
  auto printComponent = [&](const Eigen::VectorXd& eigenvalues, const Eigen::MatrixXd& eigenvectors,
                            const Eigen::Matrix3Xd& transitionDipoles) {
    auto oscillatorStrength =
        TransitionDipoleCalculator::transitionDipoleMomentToOscillatorStrength(transitionDipoles, eigenvalues);
    out << "Only principal singles CSF are printed out (CI Coef^2 > 0.1)." << std::endl;
    out << "Label is given as OccupiedOrbitalIndex -> VirtualOrbitalIndex (first orbital labelled 0)." << std::endl;
    for (int transition = 0; transition < eigenvalues.size(); ++transition) {
      double energy = eigenvalues[transition];
      const Eigen::VectorXd& state = eigenvectors.col(transition);
      const Eigen::Vector3d& trDip = transitionDipoles.col(transition);

      out << std::left << "Vertical excitation energy " << transition << ":    " << energy << " au    "
          << Utils::Constants::ev_per_hartree * energy << " eV    "
          << Utils::Constants::invCentimeter_per_hartree * energy << " cm-1    "
          << 1. / (Utils::Constants::invCentimeter_per_hartree * energy) * 1.e7 << " nm" << std::right << std::endl;

      out << std::setw(30) << "Dipole x [au]" << std::setw(30) << "Dipole y [au]" << std::setw(30) << "Dipole z [au]"
          << std::setw(30) << "Dipole Squared [au^2]" << std::setw(30) << "Oscillator Strength" << std::endl;

      out << std::setw(30) << trDip.x() << std::setw(30) << trDip.y() << std::setw(30) << trDip.z() << std::setw(30)
          << trDip.dot(trDip) << std::setw(30) << oscillatorStrength(transition) << std::endl;

      out << std::setw(30) << "MO transition labels" << std::setw(30) << "CI Coefficient" << std::setw(30)
          << "CI Coefficient squared" << std::endl;

      bool transitionsShown = false;
      for (int basisFunction = 0; basisFunction < eigenvectors.rows(); ++basisFunction) {
        if (std::pow(state(basisFunction), 2) > 0.01) {
          out << std::setw(30) << result.transitionLabels[basisFunction] << std::setw(30) << state(basisFunction)
              << std::setw(30) << state(basisFunction) * state(basisFunction) << std::endl;
          transitionsShown = true;
        }
      }
      if (!transitionsShown) {
        double maxCoef = state.maxCoeff();
        int index = 0;
        for (int i = 0; i < state.size(); ++i) {
          if (state(i) == maxCoef) {
            index = i;
            break;
          }
        }
        out << std::setw(20) << result.transitionLabels[index] << std::setw(20) << state(index) << std::setw(20)
            << state(index) * state(index) << std::endl;
      }
      out << std::endl;
    }
    out << std::endl;
  };

  if (result.singlet) {
    const Eigen::VectorXd& energies = result.singlet->eigenStates.eigenValues;
    const Eigen::MatrixXd& eigenVectors = result.singlet->eigenStates.eigenVectors;
    const Eigen::MatrixXd& transitionDipole = result.singlet->transitionDipoles;
    out << std::setw(20) << std::right << "Singlet electronic transitions." << std::endl;
    printComponent(energies, eigenVectors, transitionDipole);
  }
  if (result.triplet) {
    const Eigen::VectorXd& energies = result.triplet->eigenStates.eigenValues;
    const Eigen::MatrixXd& eigenVectors = result.triplet->eigenStates.eigenVectors;
    const Eigen::MatrixXd& transitionDipole = result.triplet->transitionDipoles;
    out << std::setw(20) << std::right << "Triplet electronic transitions." << std::endl;
    printComponent(energies, eigenVectors, transitionDipole);
  }
  if (result.unrestricted) {
    const Eigen::VectorXd& energies = result.unrestricted->eigenStates.eigenValues;
    const Eigen::MatrixXd& eigenVectors = result.unrestricted->eigenStates.eigenVectors;
    const Eigen::MatrixXd& transitionDipole = result.unrestricted->transitionDipoles;
    out << std::setw(20) << std::right << "Unrestricted electronic transitions." << std::endl;
    printComponent(energies, eigenVectors, transitionDipole);
  }
}

template<class MatrixType>
inline void matrixToCsv(std::ostream& out, const MatrixType& matrix, char delimiter) {
  const Eigen::IOFormat csvFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, std::string(1, delimiter), "\n");
  out << matrix.format(csvFormat);
  out << "\n";
}

inline Eigen::MatrixXd csvToMatrix(const std::string& filename, char delimiter) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("The file " + filename + " cannot be opened.");
  }

  std::vector<double> data;
  std::string line;
  std::string entry;
  int rowCounter = 0;

  while (std::getline(file, line)) {
    std::stringstream lineAsStream(line);
    while (std::getline(lineAsStream, entry, delimiter)) {
      data.push_back(std::stod(entry));
    }
    rowCounter++;
  }

  return Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(data.data(), rowCounter,
                                                                                            data.size() / rowCounter);
}

} // namespace Utils
} // namespace Scine

#endif // UTILS_FORMATTEDIOUTILS_H
