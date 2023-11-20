/**
 * @file
 * @brief Default module that is always loaded by Core
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "MrccState.h"
#include "Utils/IO/FilesystemHelpers.h"
#include "Utils/Technical/UniqueIdentifier.h"
#include <boost/process.hpp>

namespace Scine {
namespace Utils {
namespace ExternalQC {

MrccState::MrccState(std::string dir) : directory(std::move(dir)) {
  UniqueIdentifier id;
  stateIdentifier = id.getStringRepresentation();
  FilesystemHelpers::createDirectories(stateIdentifier);
}

MrccState::~MrccState() {
  boost::filesystem::remove_all(stateIdentifier);
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
