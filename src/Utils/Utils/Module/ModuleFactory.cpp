/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "GaussianModule.h"
#include "OrcaModule.h"

namespace Scine {
namespace Utils {

std::vector<std::shared_ptr<Scine::Core::Module>> moduleFactory() {
  return {OrcaModule::make(), GaussianModule::make()};
}

} // namespace Utils
} // namespace Scine

// At global namespace, define the entry point for the module.
#ifdef __MINGW32__
/* MinGW builds are problematic. We build with default visibility, and adding
 * an attribute __dllexport__ specifically for this singular symbol leads to the
 * loss of all other weak symbols. Essentially, here we have just expanded the
 * BOOST_DLL_ALIAS macro in order to declare the type-erased const void*
 * 'moduleFactory' without any symbol visibility attribute additions that could
 * confuse the MinGW linker, which per Boost DLL documentation is unable to mix
 * weak attributes and __dllexport__ correctly.
 *
 * If ever the default visibility for this translation unit is changed, we
 * will have to revisit this bit of code for the MinGW platform again.
 *
 * Additionally, more recent Boost releases may have fixed this problem.
 * See the macro BOOST_DLL_FORCE_ALIAS_INSTANTIATIONS as used in the library's
 * example files.
 */
extern "C" {
const void* moduleFactory = reinterpret_cast<const void*>(reinterpret_cast<intptr_t>(&Scine::Utils::moduleFactory));
}
#else
BOOST_DLL_ALIAS(Scine::Utils::moduleFactory, moduleFactory)
#endif
