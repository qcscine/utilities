/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_FORMATTEDIOUTILS_H
#define UTILS_FORMATTEDIOUTILS_H

#include <Utils/GeometricDerivatives/NormalModesContainer.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/Geometry/ElementTypes.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <iomanip>
#include <ostream>

namespace Scine {
namespace Utils {

template<class MatrixType>
std::pair<int, int> getMatrixDimensions(const MatrixType& matrix) {
  return {matrix.size(), matrix[0].size()};
};

template<>
std::pair<int, int> getMatrixDimensions(const Eigen::MatrixXd& matrix) {
  return {matrix.rows(), matrix.cols()};
};

template<>
std::pair<int, int> getMatrixDimensions(const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& matrix) {
  return {matrix.rows(), matrix.cols()};
};

template<class MatrixType>
double getElement(int row, int col, const MatrixType& matrix) {
  return matrix[row][col];
};

template<>
double getElement(int row, int col, const Eigen::MatrixXd& matrix) {
  return matrix.col(col)(row);
};

template<>
double getElement(int row, int col, const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& matrix) {
  return matrix.col(col)(row);
};

template<>
double getElement(int row, int col, const NormalModesContainer& matrix) {
  auto displCollection = matrix.getMode(col);
  return displCollection.row(row / 3)(row % 3);
};

template<class MatrixType>
void matrixPrettyPrint(std::ostream& out, const MatrixType& matrix) {
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
    out << std::endl;
  }
}

void printElement(std::ostream& out, const Eigen::SparseMatrix<double>::InnerIterator& it) {
  out << std::setw(7) << std::left << "<" + std::to_string(it.row()) + ">";
  out << std::setw(7) << std::left << "<" + std::to_string(it.col()) + ">";
  out << std::setw(7) << std::left << std::to_string(it.value());
  out << std::setw(4) << std::left << " ";
}

void matrixPrettyPrint(std::ostream& out, const Eigen::SparseMatrix<double> matrix, double threshold) {
  int nEl = 0;
  int nRows = 0;
  int headerCount = 0;
  const int nCols = 5;
  const int maxElementsPerColumn = 50;
  // print header line
  out << std::setw(10) << " ";
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
}
void matrixPrettyPrint(std::ostream& out, const NormalModesContainer& container, const ElementTypeCollection& elementTypes) {
  auto dimension = container.size();
  auto rest = dimension % 5;
  auto rows = container.getMode(0).size();
  auto columnsToRest = dimension - rest;
  auto waveNumbers = container.getWaveNumbers();

  out << std::scientific << std::right;
  // Print the matrix in blocks of 5 columns
  for (int col = 0; col < dimension - 5; col += 5) {
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

      if (row % 3 == 0)
        rowName += " x";
      else if (row % 3 == 1)
        rowName += " y";
      else
        rowName += " z";
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
      if (row % 3 == 0)
        rowName += " x";
      else if (row % 3 == 1)
        rowName += " y";
      else
        rowName += " z";
      out << std::setw(10) << rowName;
      for (int i = 0; i < rest; ++i) {
        out << std::setw(15) << std::scientific << getElement(row, columnsToRest + i, container);
      }
      out << std::endl;
    }
    out << std::defaultfloat << std::endl;
  }
}

} // namespace Utils
} // namespace Scine

#endif // UTILS_FORMATTEDIOUTILS_H
