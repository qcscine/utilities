/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_METHODINITIALIZER_H
#define UTILS_METHODINITIALIZER_H

#include <Utils/Typenames.h>
#include <string>

namespace Scine {

namespace Utils {
class AtomCollection;
} // namespace Utils

namespace Utils {

/*!
 * @brief MethodInitializer is an abstract base class for the method initializers.
 * \date 20.11.2015
 */

class MethodInitializer {
 public:
  virtual ~MethodInitializer() = default;

  /*! @brief Set the root folder where the parameters are included.
      Calling this function will have no impact on existing MethodInitializers */
  static void setRootParameterFolder(std::string path);
  /*! @brief Get the folder that contains the parameters. */
  static std::string getRootParameterFolder();

  /*! @brief Initialize the method from a Utils::AtomCollection. */
  void initializeMethod(const Utils::AtomCollection& structure);
  /*! @brief Initialize the method from Positions and atom types. */
  void initializeMethod(const Utils::PositionCollection& positions, const Utils::ElementTypeCollection& elementTypes);
  /*! @brief Initialize the method from atom types. */
  void initializeMethod(const Utils::ElementTypeCollection& elementTypes);

 protected:
  /*! initialize method to be overwritten (NVI pattern) */
  virtual void initialize(const Utils::PositionCollection& positions, const Utils::ElementTypeCollection& elementTypes) = 0;
  /*! initialize method to be overwritten (NVI pattern) */
  virtual void initialize(const Utils::ElementTypeCollection& elementTypes) = 0;

  static std::string rootParameterFolder_; //!< Default parameter path
};

} // namespace Utils
} // namespace Scine
#endif // UTILS_METHODINITIALIZER_H