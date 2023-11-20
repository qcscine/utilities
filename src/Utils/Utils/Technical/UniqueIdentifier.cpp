/**
 * @file
 * @brief A file containing definitions of classes that are just different names
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "UniqueIdentifier.h"
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

namespace Scine {
namespace Utils {

struct UniqueIdentifier::Impl {
  boost::uuids::uuid id;
};

UniqueIdentifier::~UniqueIdentifier() = default;

UniqueIdentifier::UniqueIdentifier() {
  pimpl = std::make_unique<Impl>();
  pimpl->id = boost::uuids::random_generator()();
}

UniqueIdentifier::UniqueIdentifier(const UniqueIdentifier& rhs) {
  pimpl = std::make_unique<Impl>();
  pimpl->id = rhs.pimpl->id;
}

UniqueIdentifier::UniqueIdentifier(UniqueIdentifier&& rhs) noexcept {
  pimpl = std::make_unique<Impl>();
  pimpl->id = rhs.pimpl->id;
  rhs.pimpl->id = boost::uuids::random_generator()();
}

UniqueIdentifier& UniqueIdentifier::operator=(const UniqueIdentifier& rhs) {
  pimpl->id = rhs.pimpl->id;
  return *this;
}

UniqueIdentifier& UniqueIdentifier::operator=(UniqueIdentifier&& rhs) noexcept {
  pimpl->id = rhs.pimpl->id;
  rhs.pimpl->id = boost::uuids::random_generator()();
  return *this;
}

std::string UniqueIdentifier::getStringRepresentation() const {
  return boost::uuids::to_string(pimpl->id);
}

bool UniqueIdentifier::operator==(const UniqueIdentifier& rhs) const {
  return pimpl->id == rhs.pimpl->id;
}

bool UniqueIdentifier::operator!=(const UniqueIdentifier& rhs) const {
  return !operator==(rhs);
}

bool UniqueIdentifier::operator<(const UniqueIdentifier& rhs) const {
  return pimpl->id < rhs.pimpl->id;
}

} // namespace Utils
} // namespace Scine
