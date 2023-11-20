/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Core/BaseClasses/ObjectWithLog.h>
#include <Core/Log.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine::Core;

namespace {

/* Need to make a wrapper class because otherwise pybind would need to manage
 * std::ostream smart pointers, and it's probably best to avoid that
 */
struct Sink {
  struct cout_tag {};
  struct cerr_tag {};

  Sink(const std::string& filename) : ptr(Log::fileSink(filename)) {
  }
  Sink(cout_tag /* s */) : ptr(Log::coutSink()) {
  }
  Sink(cerr_tag /* s */) : ptr(Log::cerrSink()) {
  }

  Log::SinkPtr ptr;
};

void addSink(Log::Domain& domain, std::string name, const Sink& sink) {
  domain.add(name, sink.ptr);
}

void lazy(Log::Domain& domain, std::function<std::string()> fn) {
  domain.lazy(fn);
}

} // namespace

void init_log(pybind11::module& m) {
  pybind11::class_<Sink> sink(m, "Sink", "Abstract object into which logging information can be sunk");

  pybind11::class_<Log, std::shared_ptr<Log>> logger(m, "Log");
  logger.doc() = R"(
    Multi-domain multi-sink logger

    A log object consists of multiple domains. These domains are debug, warning,
    error and output, in ascending order of importance to an end user.

    Each domain has its own set of sinks into which information can be fed.
    Multiple domains can have the same sink. Imagine a breadboard with all the
    little cables that you have to manually connect.

    By default, the debug domain has no sinks. The warning and error domains
    have a sink to stderr named "cerr". The output domain has a sink to stdout
    named "cout".

    >>> log = core.Log()
    >>> log.debug.has_sinks()
    False
    >>> log.error.has_sinks()
    True
    >>> log.debug.line("Hello world")
    >>> log.debug.add("cout", core.Log.cout_sink())
    >>> log.debug.has_sinks()
    True
    >>> log.error.remove("cerr") # Remove the stderr sink from the error domain
    >>> log.error.add("logfile", core.Log.file_sink("errors.log")) # Add a file sink instead
  )";

  logger.def(pybind11::init<>(), "Default initialize");

  logger.def_static(
      "file_sink", [](const std::string& filename) { return Sink(filename); }, "Creates a file sink");
  logger.def_static(
      "cout_sink", []() { return Sink(Sink::cout_tag{}); }, "Creates a sink to cout");
  logger.def_static(
      "cerr_sink", []() { return Sink(Sink::cerr_tag{}); }, "Creates a sink to cerr");

  logger.def_static("silent", &Log::silent, "Returns a silent log (i.e. all of its domains have no sinks)");

  pybind11::class_<Log::Domain> domain(logger, "Domain");
  domain.def(pybind11::init<>(), "Default initialize");
  domain.def("add", &addSink, "Adds a named sink to the domain");
  domain.def(
      "remove", [](Log::Domain& domain, const std::string& name) -> void { domain.remove(name); },
      "Removes a named sink from the domain.");
  domain.def("clear", &Log::Domain::clear, "Removes all sink from the domain");
  domain.def("has_sinks", &Log::Domain::operator bool);
  domain.def("line", &Log::Domain::line, "Write a line to all sinks");
  domain.def("lazy", &lazy, "Calls a string composing function and sinks its output only if the log has sinks");

  logger.def_readwrite("debug", &Log::debug, "Access the log's debug domain");
  logger.def_readwrite("warning", &Log::warning, "Access the log's warning domain");
  logger.def_readwrite("error", &Log::error, "Access the log's error domain");
  logger.def_readwrite("output", &Log::output, "Access the log's output domain");
}
