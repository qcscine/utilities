/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/IO/TurbomoleMinimalBasisfile.h"
#include "Utils/Geometry/ElementInfo.h"
#include "boost/filesystem.hpp"
#include <algorithm>
#include <boost/fusion/adapted/struct/adapt_struct.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/phoenix/fusion/at.hpp>
#include <boost/phoenix/stl/container.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/qi.hpp>
#include <iostream>

namespace Scine {
namespace Utils {
namespace detail {

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

namespace symbols {

// Mapping from element type name strings to their atomic number
struct element_type_ : qi::symbols<char, int> {
  element_type_() {
    for (const auto& nameElementPair : ElementInfo::stringToElementType()) {
      const auto& symbol = nameElementPair.first;
      /* Since we will be using these in case-insensitive matching, we need to
       * have lowercase symbols, make sure this precondition is kept
       */
      assert(std::all_of(std::begin(symbol), std::end(symbol), [](const char x) { return std::islower(x); }) &&
             "All symbols must be lowercase for case-insensitive matching!");
      if (symbol != "none") {
        add(symbol, Utils::ElementInfo::Z(nameElementPair.second));
      }
    }
  };
} const element_type;

// Mapping from spd to angular momentum number
struct angular_momentum_ : qi::symbols<char, unsigned> {
  angular_momentum_() {
    add("s", 0)("p", 1)("d", 2);
  }
} const angular_momentum;

} // namespace symbols

/* Definitions of structs that better reflect the structure of the file being
 * parsed than their actual Utils types. They're largely similar, but the Utils
 * types are a lot less easily integrable (no public members, members set by
 * functions).
 *
 * Each carries a member to convert it into the Utils type
 */
// Gtf precursor
struct GtfBase {
  double exponent;
  double coefficient;

  Gtf interpret(const int angularMomentum) const {
    return {angularMomentum, exponent, coefficient};
  }
};

// GtoExpansion precursor
struct GtoBase {
  unsigned angularMomentum;
  std::vector<GtfBase> gtfs;

  GtoExpansion interpret() const {
    GtoExpansion expansion;
    expansion.angularMomentum = angularMomentum;
    for (const GtfBase& base : gtfs) {
      expansion.gtfs.push_back(base.interpret(angularMomentum));
    }

    return expansion;
  }
};

// Map pair precursor to atomic number and AtomicGto
struct ElementAtomicGtoBase {
  int Z;
  std::vector<GtoBase> parts;

  std::pair<int, AtomicGtos> interpret() const {
    AtomicGtos atomic;

    for (const GtoBase& base : parts) {
      if (base.angularMomentum == 0) {
        atomic.s = base.interpret();
      }
      else if (base.angularMomentum == 1) {
        atomic.p = base.interpret();
      }
      else if (base.angularMomentum == 2) {
        atomic.d = base.interpret();
      }
    }

    return std::make_pair(Z, atomic);
  }
};

} // namespace detail
} // namespace Utils
} // namespace Scine

BOOST_FUSION_ADAPT_STRUCT(Scine::Utils::detail::GtfBase, (double, exponent), (double, coefficient))

BOOST_FUSION_ADAPT_STRUCT(Scine::Utils::detail::GtoBase, (unsigned, angularMomentum),
                          (std::vector<Scine::Utils::detail::GtfBase>, gtfs))

BOOST_FUSION_ADAPT_STRUCT(Scine::Utils::detail::ElementAtomicGtoBase, (int, Z),
                          (std::vector<Scine::Utils::detail::GtoBase>, parts))

namespace Scine {
namespace Utils {
namespace detail {

// Basisfile parser grammar
template<typename Iterator>
struct Basisfile : qi::grammar<Iterator> {
  // Result of parsing
  std::vector<ElementAtomicGtoBase> data;

  Basisfile() : Basisfile::base_type(start) {
    using boost::phoenix::at_c;
    using boost::phoenix::push_back;
    using qi::_1;
    using qi::_val;
    using qi::char_;
    using qi::double_;
    using qi::eol;
    using qi::lit;
    using qi::uint_;

    /* Rules */
    // clang-format off
    gtf %= qi::skip(ascii::blank)[double_ >> double_ >> eol];

    block = qi::skip(ascii::blank)[
      uint_ >> symbols::angular_momentum[at_c<0>(_val) = _1] >> eol
      >> +(gtf[push_back(at_c<1>(_val), _1)])
    ];

    parameterPair = qi::skip(ascii::blank)[
      symbols::element_type[at_c<0>(_val) = _1]
      >> +qi::char_(R"(A-Za-z0-9\-)") >> eol
      >> qi::lit("*") >> eol
      >> +(block[push_back(at_c<1>(_val), _1)])
      >> qi::lit("*") >> eol
    ];

    start = qi::skip(ascii::blank)[
      lit("$basis") >> eol
      >> lit("*") >> eol
      >> +(parameterPair[push_back(boost::phoenix::ref(data), _1)])
      >> lit("$end") >> (eol | qi::eps)
      >> qi::eoi
    ];
    // clang-format on

    /* Error handling */
    qi::on_error<qi::fail>(gtf, std::cerr << boost::phoenix::val("Expected gtf exponent and coefficient here:\"")
                                          << boost::phoenix::construct<std::string>(qi::_3, qi::_2)
                                          << boost::phoenix::val("\"\n"));

    qi::on_error<qi::fail>(block, std::cerr << boost::phoenix::val("Expected gto block here:\"")
                                            << boost::phoenix::construct<std::string>(qi::_3, qi::_2)
                                            << boost::phoenix::val("\"\n"));

    qi::on_error<qi::fail>(parameterPair, std::cerr
                                              << boost::phoenix::val("Expected element type and gto blocks here:\"")
                                              << boost::phoenix::construct<std::string>(qi::_3, qi::_2)
                                              << boost::phoenix::val("\"\n"));

    qi::on_error<qi::fail>(start, std::cerr << boost::phoenix::val("Expected basic structure of basisfile here:\"")
                                            << boost::phoenix::construct<std::string>(qi::_3, qi::_2)
                                            << boost::phoenix::val("\"\n"));
  }

  // Turn parsed precursor data into Utils types
  std::unordered_map<int, AtomicGtos> interpret() const {
    std::unordered_map<int, AtomicGtos> map;

    for (const auto& pair : data) {
      map.insert(pair.interpret());
    }

    return map;
  }

  qi::rule<Iterator, GtfBase()> gtf;
  qi::rule<Iterator, GtoBase()> block;
  qi::rule<Iterator, ElementAtomicGtoBase()> parameterPair;
  qi::rule<Iterator> start;
};

} // namespace detail

std::unordered_map<int, AtomicGtos> readTurbomoleBasisfile(const std::string& parameterFile) {
  if (!boost::filesystem::exists(parameterFile)) {
    throw std::runtime_error("File to read does not exist");
  }

  std::ifstream ifs(parameterFile);
  ifs.unsetf(std::ios::skipws);

  using Iterator = boost::spirit::istream_iterator;
  using Parser = detail::Basisfile<Iterator>;

  Iterator iter(ifs);
  Iterator end;
  Parser parser;
  bool result = boost::spirit::qi::parse(iter, end, parser);

  if (result && iter == end) {
    return parser.interpret();
  }

  throw std::runtime_error("Failed to parse basisfile");
}

} // namespace Utils
} // namespace Scine
