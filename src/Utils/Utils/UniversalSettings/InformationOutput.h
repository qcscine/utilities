/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UNIVERSALSETTINGS_INFORMATIONOUTPUT_H
#define UNIVERSALSETTINGS_INFORMATIONOUTPUT_H
/* External Headers */
#include <ostream>

namespace Scine {
namespace Utils {
namespace UniversalSettings {
/* Forward Declarations */
class DescriptorCollection;

/*!
 * This class prints information about setting descriptors in a human-readable format.
 */
class InformationOutput {
 public:
  static void print(const std::string& key, const DescriptorCollection& settings, std::ostream& out,
                    int indentation = 0, bool outputCollectionTitle = true);

  /*! Less concise-version, may contain more details than print(...). */
  static void printLong(const std::string& key, const DescriptorCollection& settings, std::ostream& out, int indentation = 0);
};

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */

#endif // UNIVERSALSETTINGS_INFORMATIONOUTPUT_H
