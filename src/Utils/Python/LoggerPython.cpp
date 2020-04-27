/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/IO/Logger.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Utils;

namespace detail {

/* These functions are all templated to admit const char* in streaming IO
 * in Log, so these helper functions instantiate that directly.
 */

void trace(const std::string& str) {
  Log::trace(str);
}

void debug(const std::string& str) {
  Log::debug(str);
}

void info(const std::string& str) {
  Log::info(str);
}

void warning(const std::string& str) {
  Log::warning(str);
}

void error(const std::string& str) {
  Log::error(str);
}

void fatal(const std::string& str) {
  Log::fatal(str);
}

} // namespace detail

void init_logger(pybind11::module& m) {
  pybind11::class_<Log> logger(m, "Log");

  pybind11::enum_<Log::severity_level> severity_level(logger, "severity_level", pybind11::arithmetic(),
                                                      "Severity level for logging");

  // clang-format off
  severity_level.value("trace", Log::severity_level::trace)
   .value("debug", Log::severity_level::debug)
   .value("info", Log::severity_level::info)
   .value("warning", Log::severity_level::warning)
   .value("error", Log::severity_level::error)
   .value("fatal", Log::severity_level::fatal);
  // clang-format on

  logger.def_static("disable_logging", &Log::disableLogging, "Disable all logging");
  logger.def_static("enable_logging", &Log::enableLogging, "Enable all logging");
  logger.def_static("start_file_logging",
                    pybind11::overload_cast<std::string, Log::severity_level, bool>(&Log::startFileLogging),
                    "Write the log to a file");
  logger.def_static("start_file_logging",
                    pybind11::overload_cast<std::string, const std::string&, bool>(&Log::startFileLogging),
                    "Write the log to a file");
  logger.def_static("start_console_logging", pybind11::overload_cast<Log::severity_level>(&Log::startConsoleLogging),
                    "Write the log to the console");
  logger.def_static("start_console_logging", pybind11::overload_cast<const std::string&>(&Log::startConsoleLogging),
                    "Write the log to the console");
  logger.def_static("stop_console_logging", &Log::stopConsoleLogging, "Stops writing the log to the console");

  logger.def_static("trace", &detail::trace, "Logs a supplied string at trace severity.");
  logger.def_static("debug", &detail::debug, "Logs a supplied string at debug severity.");
  logger.def_static("info", &detail::info, "Logs a supplied string at info severity.");
  logger.def_static("warning", &detail::warning, "Logs a supplied string at warning severity.");
  logger.def_static("error", &detail::error, "Logs a supplied string at error severity.");
  logger.def_static("fatal", &detail::fatal, "Logs a supplied string at fatal severity.");
}
