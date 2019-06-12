/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Constants.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/XYZStreamHandler.h>
#include <Utils/Technical/ScopedLocale.h>
#include <gmock/gmock.h>
#include <sstream>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {
/**
 * @class XYZStreamHandlerTest
 * @brief Comprises tests for the class Scine::Utils::AtomCollection.
 * @test
 */
class XYZStreamHandlerTest : public Test {
 public:
  std::string germanLocale;

 protected:
  void SetUp() override {
#ifdef _WIN32
    germanLocale = "de-DE";
#else
    germanLocale = "de_DE.utf8";
#endif
  }
};

TEST_F(XYZStreamHandlerTest, CorrectImportWithConversionFromAngstromToBohr) {
  std::stringstream ss("2\n\n"
                       "H      0.0  0.0  -2.0\n"
                       "V      3.0  0.4   0.0\n");

  AtomCollection structure = XYZStreamHandler::read(ss);

  auto expectedElements = ElementTypeCollection{ElementType::H, ElementType::V};
  PositionCollection expectedPositions(2, 3);
  expectedPositions << 0, 0, -2 * Constants::bohr_per_angstrom, 3 * Constants::bohr_per_angstrom,
      0.4 * Constants::bohr_per_angstrom, 0;
  ASSERT_THAT(structure.size(), Eq(2));
  ASSERT_THAT(structure.getElements(), Eq(expectedElements));
  ASSERT_THAT(structure.getPositions(), Eq(expectedPositions));
}

TEST_F(XYZStreamHandlerTest, AtomCountIsNotReliedUpon) {
  std::stringstream ss1("Dummy\n\n"
                        "H      0.0  0.0  -2.0\n"
                        "V      3.0  0.4   0.0\n");

  std::stringstream ss2("0\n\n"
                        "H      0.0  0.0  -2.0\n"
                        "V      3.0  0.4   0.0\n");

  std::stringstream ss3("10\n\n"
                        "H      0.0  0.0  -2.0\n"
                        "V      3.0  0.4   0.0\n");

  auto checkStructure = [](auto& ss) {
    auto structure = XYZStreamHandler::read(ss);
    ASSERT_THAT(structure.size(), Eq(2));
    ASSERT_THAT(structure.getPositions().rows(), Eq(2));
    ASSERT_THAT(structure.getElements().size(), Eq(2));
  };

  checkStructure(ss1);
  checkStructure(ss2);
  checkStructure(ss3);
}

TEST_F(XYZStreamHandlerTest, AllowsForSpacesAndTabsInFirstLine) {
  std::stringstream ss("  \t2 \t   \n\n"
                       "H      0.0  0.0  -2.0\n"
                       "V      3.0  0.4   0.0\n");

  ASSERT_NO_THROW(XYZStreamHandler::read(ss));
}

TEST_F(XYZStreamHandlerTest, ThrowsForInvalidLineInput) {
  std::stringstream ss1("2\n\n"
                        "H      0.0  0k0  -2.0\n"
                        "V      3.0  0.4   0.0\n");
  std::stringstream ss2("2\n\n"
                        "H      0.0  0k0  -2.0\n"
                        "V      3.0  0.4      \n");
  std::stringstream ss3("2\n\n"
                        "5      0.0  0k0  -2.0\n"
                        "V      3.0  0.4      \n");

  ASSERT_THROW(XYZStreamHandler::read(ss1), FormattedStreamHandler::FormatMismatchException);
  ASSERT_THROW(XYZStreamHandler::read(ss2), FormattedStreamHandler::FormatMismatchException);
  ASSERT_THROW(XYZStreamHandler::read(ss3), FormattedStreamHandler::FormatMismatchException);
}

TEST_F(XYZStreamHandlerTest, AllowsForFurtherColumns) {
  std::stringstream ss1("2\n\n"
                        "H      0.0  0.0  -2.0 and some text later on\n"
                        "V      3.0  0.4   0.0\n");
  std::stringstream ss2("2\n\n"
                        "H      0.0  0.0  -2.0    0.4\n"
                        "V      3.0  0.4   0.0    0.4\n");

  ASSERT_NO_THROW(XYZStreamHandler::read(ss1));
  ASSERT_NO_THROW(XYZStreamHandler::read(ss2));
}

TEST_F(XYZStreamHandlerTest, ThrowsForInvalidAtomTypes) {
  std::stringstream ss1("2\n\n"
                        "Zu     0.0  0.0  -2.0\n"
                        "Vul    3.0  0.4   0.0\n");

  ASSERT_THROW(XYZStreamHandler::read(ss1), FormattedStreamHandler::FormatMismatchException);
}

TEST_F(XYZStreamHandlerTest, AllowsLeadingSpaces) {
  // F.i. as in ORCA-generated xyz files.
  std::stringstream ss("2\n\n"
                       " H      0.0  0.0  -2.0\n"
                       "   V      3.0  0.4   0.0\n");

  auto structure = XYZStreamHandler::read(ss);

  auto expectedElements = ElementTypeCollection{ElementType::H, ElementType::V};
  PositionCollection expectedPositions(2, 3);
  expectedPositions << 0, 0, -2 * Constants::bohr_per_angstrom, 3 * Constants::bohr_per_angstrom,
      0.4 * Constants::bohr_per_angstrom, 0;

  ASSERT_THAT(structure.size(), Eq(2));
  ASSERT_THAT(structure.getElements(), Eq(expectedElements));
  ASSERT_THAT(structure.getPositions(), Eq(expectedPositions));
}

TEST_F(XYZStreamHandlerTest, AllowsForAllUppercaseSymbols) {
  // F.i. as generated in some packages
  std::stringstream ss("2\n\n"
                       "HE      0.0  0.0  -2.0\n"
                       "CA      3.0  0.4   0.0\n");

  auto structure = XYZStreamHandler::read(ss);

  auto expectedElements = ElementTypeCollection{ElementType::He, ElementType::Ca};

  ASSERT_THAT(structure.size(), Eq(2));
  ASSERT_THAT(structure.getElements(), Eq(expectedElements));
}

TEST_F(XYZStreamHandlerTest, AllowsForAllLowercaseSymbols) {
  // F.i. as generated in some packages
  std::stringstream ss("2\n\n"
                       "he      0.0  0.0  -2.0\n"
                       "ca      3.0  0.4   0.0\n");

  auto structure = XYZStreamHandler::read(ss);

  auto expectedElements = ElementTypeCollection{ElementType::He, ElementType::Ca};

  ASSERT_THAT(structure.size(), Eq(2));
  ASSERT_THAT(structure.getElements(), Eq(expectedElements));
}

TEST_F(XYZStreamHandlerTest, CorrectImportInLocaleWithCommas) {
  /* The locale de_DE.utf8 may not be present on a compiling machine. If it
   * isn't, and std::runtime_error is thrown (see std::locale constructors)
   * ignore the exception and do nothing.
   *
   * It's not possible to limit catching the runtime_error to the ScopedLocale
   * constructor since it cannot be default-constructed. Any further
   * runtime_errors caused in the remaining try-block are ignored also and
   * may hide bugs.
   */
  try {
    ScopedLocale localeWithCommas(germanLocale);
    std::stringstream ss("2\n\n"
                         "H      0.0  0.0  -2.0\n"
                         "V      3.0  0.4   0.0\n");

    auto structure = XYZStreamHandler::read(ss);

    auto expectedElements = ElementTypeCollection{ElementType::H, ElementType::V};
    PositionCollection expectedPositions(2, 3);
    expectedPositions << 0, 0, -2 * Constants::bohr_per_angstrom, 3 * Constants::bohr_per_angstrom,
        0.4 * Constants::bohr_per_angstrom, 0;

    ASSERT_THAT(structure.size(), Eq(2));
    ASSERT_THAT(structure.getElements(), Eq(expectedElements));
    ASSERT_THAT(structure.getPositions(), Eq(expectedPositions));
  }
  catch (std::runtime_error& e) {
    // Do nothing.
  }
}

std::string dimethylbutaneXYZ = "20\n"
                                "C6H14 Butane, 2,2-dimethyl- 75832\n"
                                "C          1.90460       -1.07450        0.00000\n"
                                "C          1.17020        0.26960        0.00000\n"
                                "C         -0.37650        0.22000        0.00000\n"
                                "C         -0.89830        1.66820        0.00000\n"
                                "C         -0.89830       -0.49790        1.25650\n"
                                "C         -0.89830       -0.49790       -1.25650\n"
                                "H         -1.99100        1.69140        0.00000\n"
                                "H          2.98490       -0.91340        0.00000\n"
                                "H         -0.52470       -0.01940        2.16560\n"
                                "H         -0.52470       -0.01940       -2.16560\n"
                                "H         -0.59660       -1.54710        1.28080\n"
                                "H         -0.59660       -1.54710       -1.28080\n"
                                "H         -1.99060       -0.47110        1.28940\n"
                                "H         -1.99060       -0.47110       -1.28940\n"
                                "H          1.66510       -1.67180        0.88210\n"
                                "H          1.66510       -1.67180       -0.88210\n"
                                "H          1.49110        0.84560        0.87480\n"
                                "H          1.49110        0.84560       -0.87480\n"
                                "H         -0.55200        2.21230        0.88290\n"
                                "H         -0.55200        2.21230       -0.88290";

TEST_F(XYZStreamHandlerTest, SelfConsistent) {
  std::stringstream in(dimethylbutaneXYZ), out;
  AtomCollection atoms = XYZStreamHandler::read(in);
  XYZStreamHandler::write(out, atoms);
  AtomCollection secondTime = XYZStreamHandler::read(out);

  // Compare AtomCollections
  ASSERT_THAT(atoms, Eq(secondTime));
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
