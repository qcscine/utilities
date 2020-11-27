/**
 * @file
 * @brief XYZ-formatted streaming IO
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef INCLUDE_SCINE_UTILS_IO_XYZ_STREAM_H
#define INCLUDE_SCINE_UTILS_IO_XYZ_STREAM_H

#include "Utils/IO/ChemicalFileFormats/FormattedStreamHandler.h"
#include "boost/optional/optional_fwd.hpp"
#include <fstream>

namespace Scine {

namespace Utils {

class AtomCollection;
class BondOrderCollection;

/**
 * @brief Handler for XYZ stream IO
 *
 * This class implements the FormattedStreamHandler interface for use with
 * Core's module management. The core stream IO functions are available as
 * static functions, too.
 *
 * @note Since the XYZ file format does not include bond information in any
 *   part, the read/write pair involving BondOrderCollections just throws in
 *   all cases.
 */
class XyzStreamHandler : public FormattedStreamHandler {
 public:
  //! Model string identifier for use in the Module system
  static constexpr const char* model = "XyzStreamHandler";

  //!@name Static members
  //!@{
  static AtomCollection read(std::istream& is);
  static void write(std::ostream& os, const AtomCollection& atoms);
  //!@}

  //!@name FormattedStreamHandler interface
  //!@{
  std::pair<AtomCollection, BondOrderCollection> read(std::istream& /* is */, const std::string& format) const final;

  void write(std::ostream& os, const std::string& format, const AtomCollection& atoms) const final;

  void write(std::ostream& /* os */, const std::string& format, const AtomCollection& /* atoms */,
             const BondOrderCollection& /* bondOrders */
             ) const final;

  std::vector<FormatSupportPair> formats() const final;

  bool formatSupported(const std::string& format, SupportType /* operation */) const final;

  std::string name() const final;
  //!@}
};

} // namespace Utils

} // namespace Scine

#endif
