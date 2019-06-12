/**
 * @file Logger.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_LOGGER_H_
#define UTILS_LOGGER_H_

#include <memory>
#include <sstream>
#include <string>

namespace Scine {
namespace Utils {

/**
 * @class Log Logger.h
 * @brief Static class for logging.
 *
 * By default, logging is disabled.
 * Logging can be started with startConsoleLogging() or startFileLogging().
 * It is possible to have several files for logging.
 */
class Log {
  struct ChainedLogger;

 public:
  /**
   * @brief Enum for the log entry severity.
   */
  enum class severity_level { trace, debug, info, warning, error, fatal };
  /**
   * @brief Disabled constructor: this is a static class.
   */
  Log() = delete;

  /**
   * @brief This function disables all logging.
   */
  static void disableLogging();
  /**
   * @brief This function reenables all logging if disableLogging() was called previously.
   */
  static void enableLogging();
  /**
   * @brief Flush the log messages.
   *
   * It can be that log messages are temporarily stored in a buffer. This forces output to the respective stream.
   */
  static void flush();
  /**
   * @brief Start writing log to a file.
   * @param filename The path to the file to be logged to.
   * @param minimalSeverity Log messages will be written only when their severity is at least as important as this.
   * @param autoFlush Ff true, after each log message the stream will be flushed (at higher computational cost).
   */
  static void startFileLogging(std::string filename, severity_level minimalSeverity = defaultMinimalSeverity_,
                               bool autoFlush = false);
  /**
   * @brief (Re)starts logging to the console.
   * @param minimalSeverity Log messages will be written only when their severity is at least as important as this.
   */
  static void startConsoleLogging(severity_level minimalSeverity = defaultMinimalSeverity_);
  /**
   * @brief Stops logging to the console.
   */
  static void stopConsoleLogging();

  /**
   * @brief Logging function to log statementes of the severity: trace.
   *
   * Function returning a logger. Allows direct logging with the left-shift
   * operator:
   *
   * @code Log::trace() << "Information to log"; @endcode
   *
   * @return ChainedLogger The logger to be logged with.
   */
  static ChainedLogger trace();
  /**
   * @brief Logging function to log statementes of the severity: debug.
   *
   * Function returning a logger, the function0 allow dircet logging via
   * the leftshift operator:
   * @code Log::debug() << "Information to log"; @endcode
   *
   * @return ChainedLogger The logger to be logged with.
   */
  static ChainedLogger debug();
  /**
   * @brief Logging function to log statementes of the severity: info.
   *
   * Function returning a logger, the function0 allow dircet logging via
   * the leftshift operator:
   * @code Log::info() << "Information to log"; @endcode
   *
   * @return ChainedLogger The logger to be logged with.
   */
  static ChainedLogger info();
  /**
   * @brief Logging function to log statementes of the severity: warning.
   *
   * Function returning a logger, the function0 allow dircet logging via
   * the leftshift operator:
   * @code Log::warning() << "Information to log"; @endcode
   *
   * @return ChainedLogger The logger to be logged with.
   */
  static ChainedLogger warning();
  /**
   * @brief Logging function to log statementes of the severity: error.
   *
   * Function returning a logger, the function0 allow dircet logging via
   * the leftshift operator:
   * @code Log::error() << "Information to log"; @endcode
   *
   * @return ChainedLogger The logger to be logged with.
   */
  static ChainedLogger error();
  /**
   * @brief Logging function to log statementes of the severity: fatal.
   *
   * Function returning a logger, the function0 allow dircet logging via
   * the leftshift operator:
   * @code Log::fatal() << "Information to log"; @endcode
   *
   * @return ChainedLogger The logger to be logged with.
   */
  static ChainedLogger fatal();

  /**
   * @brief Logging function to log statementes of the severity: trace.
   * @tparam String
   * @param s The statement to be logged.
   */
  template<typename String>
  static void trace(String&& s) {
    trace() << std::forward<String>(s);
  }
  /**
   * @brief Logging function to log statementes of the severity: debug.
   * @tparam String
   * @param s The statement to be logged.
   */
  template<typename String>
  static void debug(String&& s) {
    debug() << std::forward<String>(s);
  }
  /**
   * @brief Logging function to log statementes of the severity: info.
   * @tparam String
   * @param s The statement to be logged.
   */
  template<typename String>
  static void info(String&& s) {
    info() << std::forward<String>(s);
  }
  /**
   * @brief Logging function to log statementes of the severity: warning.
   * @tparam String
   * @param s The statement to be logged.
   */
  template<typename String>
  static void warning(String&& s) {
    warning() << std::forward<String>(s);
  }
  /**
   * @brief Logging function to log statementes of the severity: error.
   * @tparam String
   * @param s The statement to be logged.
   */
  template<typename String>
  static void error(String&& s) {
    error() << std::forward<String>(s);
  }
  /**
   * @brief Logging function to log statementes of the severity: fatal.
   * @tparam String
   * @param s The statement to be logged.
   */
  template<typename String>
  static void fatal(String&& s) {
    fatal() << std::forward<String>(s);
  }

 private:
  struct Impl;
  static std::unique_ptr<Impl> pImpl_;
  static const severity_level defaultMinimalSeverity_ = severity_level::info;

  /**
   * @class ChainedLogger Logger.h
   * @brief Structure allowing to chain output elements to be logged.
   *
   * The boost macro is called in the destructor to avoid creating one log entry by chained output element.
   */
  struct ChainedLogger {
    /**
     * @brief Construct a new ChainedLogger object
     * @param l The severty level of this logger.
     */
    explicit ChainedLogger(severity_level l);
    /// @brief Default move constructor.
    ChainedLogger(ChainedLogger&& t) = default;
    /// @brief Custom destructor.
    ~ChainedLogger();
    /**
     * @brief Overriding the left shift operator.
     * @tparam K The type of the message.
     * @param k The message.
     * @return std::ostream& The stream logged to.
     */
    template<typename K>
    std::ostream& operator<<(K&& k) {
      os << k;
      return os;
    }
    // @brief The severity level of this logger.
    const severity_level level;
    // @brief The stream used to log to.
    std::ostringstream os;
  };
};

} /* namespace Utils */
} /* namespace Scine */

#endif // UTILS_LOGGER_H_
