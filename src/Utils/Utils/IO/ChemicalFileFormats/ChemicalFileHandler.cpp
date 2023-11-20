/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "boost/filesystem.hpp"
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/IO/ChemicalFileFormats/MolStreamHandler.h>
#include <Utils/IO/ChemicalFileFormats/OpenBabelStreamHandler.h>
#include <Utils/IO/ChemicalFileFormats/PdbStreamHandler.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <fstream>
#include <memory>

namespace Scine {
namespace Utils {

namespace detail {

template<typename IStream, typename... Args>
auto dispatchRead(const std::string& suffix, IStream& is, Args... args) {
  std::vector<std::unique_ptr<FormattedStreamHandler>> handlers;
  handlers.emplace_back(new MolStreamHandler());
  handlers.emplace_back(new XyzStreamHandler());
  handlers.emplace_back(new PdbStreamHandler());
  handlers.emplace_back(new OpenBabelStreamHandler());

  for (auto& handlerPtr : handlers) {
    if (handlerPtr->formatSupported(suffix, FormattedStreamHandler::SupportType::ReadOnly)) {
      return handlerPtr->read(is, suffix, args...);
    }
  }

  // If no handler is found that matches the suffix, it is unsupported
  throw FormattedStreamHandler::FormatUnsupportedException();
}

template<typename OStream, typename... Args>
auto dispatchWrite(const std::string& suffix, OStream& os, Args... args) {
  std::vector<std::unique_ptr<FormattedStreamHandler>> handlers;
  handlers.emplace_back(new MolStreamHandler());
  handlers.emplace_back(new XyzStreamHandler());
  handlers.emplace_back(new PdbStreamHandler());
  handlers.emplace_back(new OpenBabelStreamHandler());

  for (auto& handlerPtr : handlers) {
    if (handlerPtr->formatSupported(suffix, FormattedStreamHandler::SupportType::WriteOnly)) {
      handlerPtr->write(os, suffix, args...);
      return;
    }
  }

  // If no handler is found that matches the suffix, it is unsupported
  throw FormattedStreamHandler::FormatUnsupportedException();
}

std::string getSuffix(const boost::filesystem::path& filepath) {
  // Match the file suffix against supported formats
  std::string suffix = filepath.extension().string();

  // The suffix may be absent or just a dot, hence we cannot deduce a format
  if (suffix.size() <= 1) {
    throw FormattedStreamHandler::FormatUnsupportedException();
  }

  // Remove the leading dot
  return suffix.substr(1);
}

} // namespace detail

std::pair<AtomCollection, BondOrderCollection> ChemicalFileHandler::read(const std::string& filename) {
  boost::filesystem::path filepath(filename);

  if (!boost::filesystem::exists(filepath)) {
    throw FileInaccessibleException();
  }

  // Try to open the file
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw FileInaccessibleException();
  }

  auto data = detail::dispatchRead(detail::getSuffix(filepath), file);

  file.close();

  return data;
}

void ChemicalFileHandler::write(const std::string& filename, const AtomCollection& atoms, const std::string& comment) {
  boost::filesystem::path filepath(filename);

  std::ofstream file(filename);
  if (!file.is_open()) {
    throw FileInaccessibleException();
  }

  detail::dispatchWrite(detail::getSuffix(filepath), file, atoms, comment);

  file.close();
}

void ChemicalFileHandler::write(const std::string& filename, const AtomCollection& atoms,
                                const BondOrderCollection& bondOrders, const std::string& comment) {
  boost::filesystem::path filepath(filename);

  std::ofstream file(filename);
  if (!file.is_open()) {
    throw FileInaccessibleException();
  }

  detail::dispatchWrite(detail::getSuffix(filepath), file, atoms, bondOrders, comment);

  file.close();
}

} // namespace Utils
} // namespace Scine
