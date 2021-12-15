/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_CALCULATION_ROUTINES_H_
#define UTILS_CALCULATION_ROUTINES_H_

#include <Core/Exceptions.h>
#include <Core/Interfaces/Calculator.h>
#include <Core/Log.h>
#include <Utils/CalculatorBasics/Results.h>
#include <boost/exception/diagnostic_information.hpp>
#include <string>

namespace Scine {
namespace Utils {
namespace CalculationRoutines {

/**
 * @brief Call calculator's calculate function with try catch and throw specific error for failure
 *
 * @param calculator A Scine calculator
 * @param log A Scine logger
 * @param errorMessage The error message of the exception
 * @throw Core::UnsuccessfulCalculationException if calculation fails
 * @return Results A Scine Results object
 */
inline Results calculateWithCatch(Core::Calculator& calculator, Core::Log& log, const std::string& errorMessage) {
  Utils::Results results;
  try {
    results = calculator.calculate("");
    if (!results.get<Property::SuccessfulCalculation>()) {
      throw Core::UnsuccessfulCalculationException("Calculator signalled unsuccessful calculation.");
    }
  }
  catch (...) {
    log.error << errorMessage << Core::Log::nl;
    throw;
  }
  return results;
}

/**
 *  * @brief Give calculator a new logger based on given booleans
 */
inline void setLog(Core::Calculator& calculator, bool error, bool warning, bool out) {
  auto log = Core::Log::silent();
  if (out) {
    log.output.add("cout", Core::Log::coutSink());
  }
  if (warning) {
    log.warning.add("cerr", Core::Log::cerrSink());
  }
  if (error) {
    log.error.add("cerr", Core::Log::cerrSink());
  }
  calculator.setLog(log);
}

} // namespace CalculationRoutines
} // namespace Utils
} // namespace Scine

#endif // UTILS_CALCULATION_ROUTINES_H_
