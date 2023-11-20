/**
 * @file
 * @brief A file containing definitions of classes that are just different names
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_TYPENAMES_H_
#define UTILS_TYPENAMES_H_

#include "Utils/Geometry/ElementTypes.h"
#include <Eigen/Core>
#include <cassert>
#include <vector>

namespace Scine {
namespace Utils {

/**
 * @class Scine::Utils::Displacement Typenames.h
 * @brief Another name for an Eigen::RowVector3d .
 */
using Displacement = Eigen::RowVector3d;
/**
 * @class Scine::Utils::DisplacementCollection Typenames.h
 * @brief Another name for a row-major Eigen::MatrixX3d .
 */
using DisplacementCollection = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;

/**
 * @class Scine::Utils::ElementTypeCollection Typenames.h
 * @brief Another name for an std::vector<ElementType> .
 */
using ElementTypeCollection = std::vector<ElementType>;

/**
 * @class Scine::Utils::ElementTypeCollection Typenames.h
 * @brief Another name for an std::pair<std::tuple<std::string, std::string, unsigned int>>.
 */
using ResidueInformation = std::tuple<std::string, std::string, unsigned int>;

/**
 * @class Scine::Utils::ElementTypeCollection Typenames.h
 * @brief Another name for an std::vector<ResidueInformation>.
 */
using ResidueCollection = std::vector<ResidueInformation>;

/**
 * @class Scine::Utils::Gradient Typenames.h
 * @brief Another name for an Eigen::RowVector3d .
 */
using Gradient = Eigen::RowVector3d;

/**
 * @class Scine::Utils::GradientCollection Typenames.h
 * @brief Another name for a row-major Eigen::MatrixX3d.
 */
using GradientCollection = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;

/**
 * @class Scine::Utils::HessianMatrix Typenames.h
 * @brief Another name for a Eigen::MatrixXd .
 */
using HessianMatrix = Eigen::MatrixXd;

/**
 * @class Scine::Utils::Dipole Typenames.h
 * @brief Another name for an Eigen::RowVector3d .
 */
using Dipole = Eigen::RowVector3d;

/**
 * @class Scine::Utils::DipoleGradient Typenames.h
 * @brief Another name for a row-major Eigen::MatrixX3d.
 */
using DipoleGradient = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;

/**
 * @class Scine::Utils::Position Typenames.h
 * @brief Another name for an Eigen::RowVector3d .
 */
using Position = Eigen::RowVector3d;

/**
 * @class Scine::Utils::PositionCollection Typenames.h
 * @brief Another name for a row-major Eigen::MatrixX3d .
 */
using PositionCollection = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;

} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_TYPENAMES_H_
