/**
 * @file
 * @brief Implements a StreamHandler leveraging OpenBabel for access to more
 *   file formats
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SCINE_OPENBABEL_FILEHANDLER_H
#define SCINE_OPENBABEL_FILEHANDLER_H

/* External Includes */
#include "Utils/IO/ChemicalFileFormats/FormattedStreamHandler.h"
#include <string>

namespace Scine {
namespace Utils {

/**
 * @brief If obabel is in the PATH, supports IO to range of OpenBabel's file
 *   formats
 *
 * @note Concerning the GPLv2 OpenBabel license, this StreamHandler very much
 *   communicates 'at arms length' with OpenBabel
 *   (https://www.gnu.org/licenses/old-licenses/gpl-2.0-faq.en.html#TOCGPLInProprietarySystem)
 *   by merely calling it via system process-level interaction, if present.
 *
 *   The StreamHandler has been tested to work with the following versions of OpenBabel: 2.4.1,
 *   3.1.0, 3.1.1. It does not work with OpenBabel 2.4.0.
 */
class OpenBabelStreamHandler : public FormattedStreamHandler {
 public:
  static constexpr const char* model = "OpenBabelStreamHandler";

  /**
   * @brief List of checked and supported file formats
   * @note This list is po
   */
  static const std::vector<FormatSupportPair>& getSupportedFormats();

  /**
   * @brief Checks for the presence of 'obabel' in the PATH
   *
   * @returns Whether 'obabel' is found and can be used
   */
  static bool checkForBinary();

  /**
   * @brief Converts between two formats using openbabel with two streams
   *
   * @param is The input stream for obabel
   * @param os The output stream for obabel
   * @param inFormat The format string of @p is
   * @param outFormat The format string of @p os
   *
   * @note For supported formats, @see openBabelFormats
   *
   * @warning This does not check whether obabel is present. Do not use this
   *   function without a pre-pended @p checkForBinary call.
   *
   * @pre @p is populated with data of format @p inFormat
   * @post @p os is populated with data of format @p outFormat
   *
   * @returns The exit code of obabel execution
   */
  static int indirect(std::istream& is, std::ostream& os, const std::string& inFormat, const std::string& outFormat);

  std::pair<AtomCollection, BondOrderCollection> read(std::istream& is, const std::string& format) const final;

  void write(std::ostream& os, const std::string& format, const AtomCollection& atoms, const std::string& comment = "") const final;

  void write(std::ostream& os, const std::string& format, const AtomCollection& atoms,
             const BondOrderCollection& bondOrders, const std::string& comment = "") const final;

  std::vector<FormatSupportPair> formats() const final;

  bool formatSupported(const std::string& format, SupportType operation) const final;

  std::string name() const final;

 private:
  bool _enabled = checkForBinary();
};

} // namespace Utils
} // namespace Scine

#endif /* DFT_DFTCALCULATOR_H_ */
