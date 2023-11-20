/**
 * @file
 * @brief Default module that is always loaded by Core
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILSOS_MRCCSTATE_H
#define UTILSOS_MRCCSTATE_H

#include <Core/BaseClasses/StateHandableObject.h>
#include <string>

namespace Scine::Utils::ExternalQC {

/**
 * @class
 * @brief The MRCC state.
 */
struct MrccState final : public Core::State {
  MrccState(std::string dir);
  /// @brief Destructor, deletes all files.
  ~MrccState() final;

  std::string directory;
  std::string stateIdentifier;
};

} // namespace Scine::Utils::ExternalQC

#endif // UTILSOS_MRCCSTATE_H
