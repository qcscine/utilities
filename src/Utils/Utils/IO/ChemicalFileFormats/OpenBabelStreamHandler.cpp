/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "OpenBabelStreamHandler.h"
#include "Utils/Bonds/BondOrderCollection.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/IO/ChemicalFileFormats/MolStreamHandler.h"
#include "Utils/IO/ChemicalFileFormats/XyzStreamHandler.h"
#include "boost/optional.hpp"
#include "boost/process/child.hpp"
#include "boost/process/io.hpp"
#include "boost/process/search_path.hpp"

namespace Scine {

namespace Utils {

constexpr const char* OpenBabelStreamHandler::model;

const std::vector<OpenBabelStreamHandler::FormatSupportPair>& OpenBabelStreamHandler::getSupportedFormats() {
  static std::vector<FormatSupportPair> formats;

  // If obabel isn't present or formats are already determined, exit early
  if (!checkForBinary() || !formats.empty()) {
    return formats;
  }

  /* Now we are in construct on first-use land. Call obabel to reveal the
   * supported read and write formats:
   */
  const std::string readFormatsCallString = "obabel -L formats read";
  const std::string writeFormatsCallString = "obabel -L formats write";

  boost::process::pstream readFormatsStream;
  boost::process::pstream writeFormatsStream;

  boost::process::child readFormatsChild(readFormatsCallString, boost::process::std_out > readFormatsStream);
  boost::process::child writeFormatsChild(writeFormatsCallString, boost::process::std_out > writeFormatsStream);

  readFormatsChild.wait();
  writeFormatsChild.wait();

  // Function to parse the format string identifiers from each list
  auto parseFormats = [](std::istream& is) -> std::vector<std::string> {
    std::vector<std::string> formats;

    std::string formatName;
    while (!is.eof()) {
      is >> formatName;
      if (is.fail() && !is.eof()) {
        throw FormattedStreamHandler::FormatMismatchException();
      }
      else if (is.fail() && is.eof()) {
        break;
      }

      formats.push_back(std::move(formatName));

      // Ignore the rest of the line
      is.ignore(std::numeric_limits<std::streamsize>::max(), is.widen('\n'));
    }

    return formats;
  };

  std::vector<std::string> readFormats = parseFormats(readFormatsStream);
  std::vector<std::string> writeFormats = parseFormats(writeFormatsStream);

  std::sort(std::begin(readFormats), std::end(readFormats));

  std::sort(std::begin(writeFormats), std::end(writeFormats));

  formats.reserve(readFormats.size() + writeFormats.size());

  /* Step through both lists simultaneously to determine the venn diagram of
   * only in the read list / in both lists / only in the write list.
   */
  auto readIter = std::begin(readFormats);
  const auto readEnd = std::end(readFormats);
  auto writeIter = std::begin(writeFormats);
  const auto writeEnd = std::end(writeFormats);

  while (readIter != readEnd && writeIter != writeEnd) {
    if (*readIter == *writeIter) {
      formats.emplace_back(std::move(*readIter), SupportType::ReadWrite);
      ++readIter;
      ++writeIter;
    }
    else if (*readIter < *writeIter) {
      formats.emplace_back(std::move(*readIter), SupportType::ReadOnly);
      ++readIter;
    }
    else {
      // Remaining case is *readIter > *writeIter
      formats.emplace_back(std::move(*writeIter), SupportType::WriteOnly);
      ++writeIter;
    }
  }

  /* Once either iterator hits the end, the loop above is terminated, but
   * there may be entries in read or write still
   */
  if (readIter == readEnd && writeIter != writeEnd) {
    while (writeIter != writeEnd) {
      formats.emplace_back(std::move(*writeIter), SupportType::WriteOnly);
      ++writeIter;
    }
  }

  if (readIter != readEnd && writeIter == writeEnd) {
    while (readIter != readEnd) {
      formats.emplace_back(std::move(*readIter), SupportType::ReadOnly);
      ++readIter;
    }
  }

  // We reserved too much space to avoid reallocations, so we can shrink it now
  formats.shrink_to_fit();
  return formats;
};

bool OpenBabelStreamHandler::checkForBinary() {
  return !boost::process::search_path("obabel").empty();
}

int OpenBabelStreamHandler::indirect(std::istream& is, std::ostream& os, const std::string& inFormat,
                                     const std::string& outFormat) {
  // Generate the call string
  std::string callString = "obabel -i" + inFormat + " -o" + outFormat;

  // If the input format is SMILES or InChI, add some parameters
  if (inFormat == "smi" || inFormat == "smiles" || inFormat == "inchi" || inFormat == "can") {
    // Generate three dimensional coordinates and make all hydrogens explicit
    callString += " --gen3D -h";
  }

  // Construct pipe streams for redirection
  boost::process::opstream ips;
  boost::process::pstream ps;
  boost::process::pstream err;

  // Start the child process
  boost::process::child childProcess(callString, boost::process::std_in<ips, boost::process::std_out> ps,
                                     boost::process::std_err > err);

  // Feed our istream into the child process' stdin
  ips << is.rdbuf();
  ips.flush();
  ips.pipe().close();

  // Wait for the child process to exit
  childProcess.wait();

  std::stringstream stderrStream;
#if BOOST_VERSION >= 107000
  /* NOTE: This implementation of buffer transfers in boost process has a bug
   * that isn't fixed before Boost 1.70.
   */
  os << ps.rdbuf();
  stderrStream << err.rdbuf();
#else
  // Workaround: cast to a parent class implementing rdbuf() correctly.
  using BasicIOSReference = std::basic_ios<char, std::char_traits<char>>&;
  // Feed the results into our ostream
  os << static_cast<BasicIOSReference>(ps).rdbuf();
  stderrStream << static_cast<BasicIOSReference>(err).rdbuf();
#endif

  /* Interpret the error string, it must contain (besides any warnings) at least
   * the string '1 molecule converted' in order to possibly be correct
   */
  std::string errorString = stderrStream.str();
  if (errorString.find("1 molecule converted") == std::string::npos) {
    return 1;
  }

  /* Very little reason to look at exit codes as obabel always yields 0, but
   * who knows, maybe sometime in the future?
   */
  return childProcess.exit_code();
}

std::pair<AtomCollection, BondOrderCollection> OpenBabelStreamHandler::read(std::istream& is, const std::string& format) const {
  if (!_enabled || !formatSupported(format, SupportType::ReadOnly)) {
    throw FormattedStreamHandler::FormatUnsupportedException();
  }

  /* Here we need to use the MOL format reader/writer to carry bond order
   * information
   */
  std::stringstream intermediate;

  // Call openbabel to convert is into mol format in intermediate
  int returnValue = indirect(is, intermediate, format, "mol");
  if (returnValue != 0) {
    throw FormattedStreamHandler::FormatMismatchException();
  }

  return MolStreamHandler::read(intermediate);
}

void OpenBabelStreamHandler::write(std::ostream& os, const std::string& format, const AtomCollection& atoms) const {
  if (!_enabled || !formatSupported(format, SupportType::WriteOnly)) {
    throw FormattedStreamHandler::FormatUnsupportedException();
  }

  /* Pipe AtomCollections using the xyz reader/writer to reduce extraneous
   * character parsing / writing
   */
  std::stringstream intermediate;
  XyzStreamHandler::write(intermediate, atoms);
  intermediate << EOF;

  int returnValue = indirect(intermediate, os, "xyz", format);
  if (returnValue != 0) {
    throw FormattedStreamHandler::FormatMismatchException();
  }
}

void OpenBabelStreamHandler::write(std::ostream& os, const std::string& format, const AtomCollection& atoms,
                                   const BondOrderCollection& bondOrders) const {
  if (!_enabled || !formatSupported(format, SupportType::WriteOnly)) {
    throw FormattedStreamHandler::FormatUnsupportedException();
  }

  std::stringstream intermediate;
  /* Here we need to use the MOL format reader/writer to carry bond order
   * information
   */
  MolStreamHandler::write(intermediate, atoms, bondOrders, "V2000");
  intermediate << EOF;

  int returnValue = indirect(intermediate, os, "mol", format);
  if (returnValue != 0) {
    throw FormattedStreamHandler::FormatMismatchException();
  }
}

std::vector<OpenBabelStreamHandler::FormatSupportPair> OpenBabelStreamHandler::formats() const {
  if (_enabled) {
    return getSupportedFormats();
  }

  return {};
}

bool OpenBabelStreamHandler::formatSupported(const std::string& format, SupportType process) const {
  const auto& formats = getSupportedFormats();
  auto findIter = std::find_if(std::begin(formats), std::end(formats),
                               [&format](const auto& pair) -> bool { return pair.first == format; });

  // If the format string isn't found, the format isn't supported
  if (findIter == std::end(formats)) {
    return false;
  }

  // If the format is read-write, it always is
  if (findIter->second == SupportType::ReadWrite) {
    return true;
  }

  // Otherwise, the process has to match the SupportType
  return findIter->second == process;
}

std::string OpenBabelStreamHandler::name() const {
  return OpenBabelStreamHandler::model;
}

} // namespace Utils

} // namespace Scine
