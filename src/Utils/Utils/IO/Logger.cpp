/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/IO/Logger.h"
#include <Utils/UniversalSettings/Exceptions.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <boost/core/null_deleter.hpp>
#include <boost/log/core/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/expressions/formatters/date_time.hpp>
#include <boost/log/sinks/sync_frontend.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <fstream>
#include <ostream>
#include <utility>

namespace logging = boost::log;
namespace expr = boost::log::expressions;
namespace sinks = boost::log::sinks;
namespace attrs = boost::log::attributes;

namespace Scine {
namespace Utils {

BOOST_LOG_ATTRIBUTE_KEYWORD(line_id, "LineID", unsigned int)
BOOST_LOG_ATTRIBUTE_KEYWORD(timestamp, "TimeStamp", boost::posix_time::ptime)
BOOST_LOG_ATTRIBUTE_KEYWORD(severity, "Severity", logging::trivial::severity_level)

using text_sink = boost::log::sinks::synchronous_sink<boost::log::sinks::text_ostream_backend>;
using text_sink_ptr = boost::shared_ptr<text_sink>;

namespace {
boost::log::trivial::severity_level getBoostSeverity(Log::severity_level s) {
  return static_cast<boost::log::trivial::severity_level>(s);
}
} // namespace

struct Log::Impl {
  void disableLogging();
  void enableLogging();
  void flush();
  void startFileLogging(std::string filename, severity_level minimalSeverity, bool autoFlush = false);
  void startConsoleLogging(severity_level minimalSeverity);
  void stopConsoleLogging();
  text_sink_ptr addSink(const boost::shared_ptr<std::ostream>& logStream, severity_level minimalSeverity,
                        bool autoFlush = false);
  void initialize();
  bool isActive() const {
    return enabled_ && (fileLoggingActive_ || consoleLoggingActive_);
  }
  boost::log::sources::severity_logger_mt<boost::log::trivial::severity_level> logger_;
  text_sink_ptr consoleSink_;
  bool enabled_ = true;
  bool initialized_ = false;
  bool fileLoggingActive_ = false;
  bool consoleLoggingActive_ = false;
};
std::unique_ptr<Log::Impl> Log::pImpl_ = std::make_unique<Log::Impl>();
std::mutex Log::m_;

void Log::Impl::disableLogging() {
  enabled_ = false;
}

void Log::Impl::enableLogging() {
  enabled_ = true;
}

void Log::Impl::flush() {
  logging::core::get()->flush();
}

void Log::Impl::startFileLogging(std::string filename, severity_level minimalSeverity, bool autoFlush) {
  if (!initialized_)
    initialize();
  fileLoggingActive_ = true;

  // add a logfile stream to our sink
  auto fileStream = boost::make_shared<std::ofstream>(filename);
  addSink(fileStream, minimalSeverity, autoFlush);
}

void Log::Impl::startConsoleLogging(severity_level minimalSeverity) {
  if (!initialized_)
    initialize();
  stopConsoleLogging();

  // add "console" output stream to our sink
  auto consoleStream = boost::shared_ptr<std::ostream>(&std::clog, boost::null_deleter());
  consoleSink_ = addSink(consoleStream, minimalSeverity);
  consoleLoggingActive_ = true;
}

void Log::Impl::stopConsoleLogging() {
  consoleLoggingActive_ = false;
  if (consoleSink_) {
    auto core = logging::core::get();
    // Remove the sink from the core, so that no records are passed to it
    core->remove_sink(consoleSink_);

    // Break the feeding loop
    // consoleSink_->stop(); // DOES NOT EXIST

    // Flush all log records that may have left buffered
    consoleSink_->flush();
    consoleSink_.reset();
  }
}

text_sink_ptr Log::Impl::addSink(const boost::shared_ptr<std::ostream>& logStream, severity_level minimalSeverity,
                                 bool autoFlush) {
  auto core = logging::core::get();

  // add a text sink
  using text_sink = sinks::synchronous_sink<sinks::text_ostream_backend>;
  boost::shared_ptr<text_sink> sink = boost::make_shared<text_sink>();

  sink->locked_backend()->add_stream(logStream);
  sink->locked_backend()->auto_flush(autoFlush);

  // specify the format of the log message
  logging::formatter formatter = expr::stream << std::setw(7) << std::setfill('0') << line_id << std::setfill(' ')
                                              << " | " << expr::format_date_time(timestamp, "%Y-%m-%d, %H:%M:%S.%f") << " "
                                              << "[" << logging::trivial::severity << "]"
                                              << " - " << expr::smessage;
  sink->set_formatter(formatter);

  // only messages with severity >= minimalSeverity are written
  sink->set_filter(severity >= getBoostSeverity(minimalSeverity));

  // "register" our sink
  core->add_sink(sink);
  return sink;
}

void Log::Impl::initialize() {
  initialized_ = false;
  logger_.add_attribute("LineID", attrs::counter<unsigned int>(1)); // lines are sequentially numbered
  logger_.add_attribute("TimeStamp", attrs::local_clock());         // each log line gets a timestamp
}

void Log::disableLogging() {
  std::lock_guard<std::mutex> lock(m_);
  pImpl_->disableLogging();
}

void Log::enableLogging() {
  std::lock_guard<std::mutex> lock(m_);
  pImpl_->enableLogging();
}

void Log::flush() {
  std::lock_guard<std::mutex> lock(m_);
  pImpl_->flush();
}

void Log::startFileLogging(std::string filename, severity_level minimalSeverity, bool autoFlush) {
  std::lock_guard<std::mutex> lock(m_);
  pImpl_->startFileLogging(std::move(filename), minimalSeverity, autoFlush);
}

void Log::startFileLogging(std::string filename, const std::string& minimalSeverity, bool autoFlush) {
  std::lock_guard<std::mutex> lock(m_);
  minimalSeverity.empty() ? pImpl_->startFileLogging(std::move(filename), defaultMinimalSeverity_)
                          : pImpl_->startFileLogging(std::move(filename), stringToSeverity(minimalSeverity), autoFlush);
}

void Log::startConsoleLogging(const std::string& minimalSeverity) {
  if (minimalSeverity == SettingsNames::LogLevels::none) {
    Utils::Log::stopConsoleLogging();
    return;
  }
  else {
    startConsoleLogging(stringToSeverity(minimalSeverity));
  }
}
Log::severity_level Log::stringToSeverity(const std::string& minimalSeverity) {
  severity_level requiredSeverity;
  if (minimalSeverity == SettingsNames::LogLevels::trace)
    requiredSeverity = severity_level::trace;
  else if (minimalSeverity == SettingsNames::LogLevels::debug)
    requiredSeverity = severity_level::debug;
  else if (minimalSeverity == SettingsNames::LogLevels::info)
    requiredSeverity = severity_level::info;
  else if (minimalSeverity == SettingsNames::LogLevels::warning)
    requiredSeverity = severity_level::warning;
  else if (minimalSeverity == SettingsNames::LogLevels::error)
    requiredSeverity = severity_level::error;
  else if (minimalSeverity == SettingsNames::LogLevels::fatal)
    requiredSeverity = severity_level::fatal;
  else
    throw UniversalSettings::OptionDoesNotExistException(minimalSeverity, Utils::SettingsNames::loggerVerbosity);
  return requiredSeverity;
}
void Log::startConsoleLogging(severity_level minimalSeverity) {
  std::lock_guard<std::mutex> lock(m_);
  pImpl_->startConsoleLogging(minimalSeverity);
}

void Log::stopConsoleLogging() {
  std::lock_guard<std::mutex> lock(m_);
  pImpl_->stopConsoleLogging();
}

Log::ChainedLogger Log::trace() {
  return ChainedLogger(severity_level::trace);
}
Log::ChainedLogger Log::debug() {
  return ChainedLogger(severity_level::debug);
}
Log::ChainedLogger Log::info() {
  return ChainedLogger(severity_level::info);
}
Log::ChainedLogger Log::warning() {
  return ChainedLogger(severity_level::warning);
}
Log::ChainedLogger Log::error() {
  return ChainedLogger(severity_level::error);
}
Log::ChainedLogger Log::fatal() {
  return ChainedLogger(severity_level::fatal);
}

Log::ChainedLogger::ChainedLogger(severity_level l) : level(l) {
}

Log::ChainedLogger::~ChainedLogger() {
  if (pImpl_->isActive())
    BOOST_LOG_SEV(pImpl_->logger_, getBoostSeverity(level)) << os.str();
}

} /* namespace Utils */
} /* namespace Scine */
