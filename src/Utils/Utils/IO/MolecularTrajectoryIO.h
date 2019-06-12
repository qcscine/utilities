/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_MOLECULARTRAJECTORYIO_H
#define UTILS_MOLECULARTRAJECTORYIO_H

#include "Utils/Typenames.h"
#include <string>

namespace Scine {
namespace Utils {
class MolecularTrajectory;
enum class ElementType;

/**
 * @class MolecularTrajectoryIO MolecularTrajectoryIO.h
 * @brief Class for input and output of MolecularTrajectory classes.
 */
class MolecularTrajectoryIO {
 public:
  /**
   * @brief Enum class for possible input and output formats.
   */
  enum class format { xyz, binary };
  /**
   * @brief Write a molecular trajectory to a file.
   *
   * @param f The format to be used/expected.
   * @param fileName The file path.
   * @param m The trajectory.
   * @throws std::runtime_error If the file could not be created.
   */
  static void write(format f, const std::string& fileName, const MolecularTrajectory& m);
  /**
   * @brief Write a molecular trajectory to a stream.
   *
   * @param f The format to be used/expected.
   * @param out  The output stream.
   * @param m The trajectory.
   */
  static void write(format f, std::ostream& out, const MolecularTrajectory& m);
  /**
   * @brief Read a molecular trajectory from a file.
   *
   * @param f The format to be used/expected.
   * @param fileName The file path.
   * @throws std::runtime_error If the file could not be opened.
   * @return MolecularTrajectory Returns the trajectory.
   */
  static MolecularTrajectory read(format f, const std::string& fileName);
  /**
   * @brief Read a molecular trajectory from a stream.
   *
   * @param f The format to be used/expected.
   * @param in The input stream.
   * @throws std::runtime_error If the format isn't supported.
   * @return MolecularTrajectory Returns the trajectory.
   */
  static MolecularTrajectory read(format f, std::istream& in);
  /**
   * @brief Output a single line (one atom) of a XYZ file.
   * @param out The output stream to be written to.
   * @param e The ElementType of the atom tho be printed.
   * @param p The Position of the atom tho be printed.
   */
  static void writeXYZLine(std::ostream& out, ElementType e, const Position& p);

 private:
  static void writeBinary(std::ostream& out, const MolecularTrajectory& m);
  static void writeXYZ(std::ostream& out, const MolecularTrajectory& m);
  static MolecularTrajectory readBinary(std::istream& in);
  static MolecularTrajectory readXYZ(std::istream& in);
};

} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_MOLECULARTRAJECTORYIO_H
