/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Constants.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Technical/ScopedLocale.h>
#include <gmock/gmock.h>
#include <boost/dll.hpp>
#include <boost/filesystem.hpp>
#include <sstream>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {
/**
 * @class XyzStreamHandlerTest
 * @brief Comprises tests for the class Scine::Utils::AtomCollection.
 * @test
 */
class XyzStreamHandlerTest : public Test {
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

TEST_F(XyzStreamHandlerTest, CorrectImportWithConversionFromAngstromToBohr) {
  std::stringstream ss("2\n\n"
                       "H      0.0  0.0  -2.0\n"
                       "V      3.0  0.4   0.0\n");

  AtomCollection structure = XyzStreamHandler::read(ss);

  auto expectedElements = ElementTypeCollection{ElementType::H, ElementType::V};
  PositionCollection expectedPositions(2, 3);
  expectedPositions << 0, 0, -2 * Constants::bohr_per_angstrom, 3 * Constants::bohr_per_angstrom,
      0.4 * Constants::bohr_per_angstrom, 0;
  ASSERT_THAT(structure.size(), Eq(2));
  ASSERT_THAT(structure.getElements(), Eq(expectedElements));
  ASSERT_THAT(structure.getPositions(), Eq(expectedPositions));
}

TEST_F(XyzStreamHandlerTest, ThrowsForWrongAtomCount) {
  std::stringstream ss1("1\n\n"
                        "H      0.0  0.0  -2.0\n"
                        "V      3.0  0.4   0.0\n");

  std::stringstream ss2("10\n\n"
                        "H      0.0  0.0  -2.0\n"
                        "V      3.0  0.4   0.0\n");

  ASSERT_THROW(XyzStreamHandler::read(ss1), FormattedStreamHandler::AtomNumberMismatchException);
  ASSERT_THROW(XyzStreamHandler::read(ss2), FormattedStreamHandler::AtomNumberMismatchException);
}

TEST_F(XyzStreamHandlerTest, ThrowsForWrongHeader) {
  std::stringstream ss1("Dummy\n\n"
                        "H      0.0  0.0  -2.0\n"
                        "V      3.0  0.4   0.0\n");

  std::stringstream ss2(""
                        "H      0.0  0.0  -2.0\n"
                        "V      3.0  0.4   0.0\n");

  std::stringstream ss3(" \n\n"
                        "H      0.0  0.0  -2.0\n"
                        "V      3.0  0.4   0.0\n");

  std::stringstream ss4("12.1\n\n"
                        "H      0.0  0.0  -2.0\n"
                        "V      3.0  0.4   0.0\n");

  std::stringstream ss5("-2\n\n"
                        "H      0.0  0.0  -2.0\n"
                        "V      3.0  0.4   0.0\n");

  ASSERT_THROW(XyzStreamHandler::read(ss1), FormattedStreamHandler::FormatMismatchException);
  ASSERT_THROW(XyzStreamHandler::read(ss2), FormattedStreamHandler::FormatMismatchException);
  ASSERT_THROW(XyzStreamHandler::read(ss3), FormattedStreamHandler::FormatMismatchException);
  ASSERT_THROW(XyzStreamHandler::read(ss4), FormattedStreamHandler::FormatMismatchException);
  ASSERT_THROW(XyzStreamHandler::read(ss5), FormattedStreamHandler::FormatMismatchException);
}

TEST_F(XyzStreamHandlerTest, ThrowsForStringInCoordinate) {
  std::stringstream ss1("2\n\n"
                        "H      0.0  0.0  -2.0\n"
                        "V      3.0  something   0.0\n");

  ASSERT_THROW(XyzStreamHandler::read(ss1), FormattedStreamHandler::FormatMismatchException);
}

TEST_F(XyzStreamHandlerTest, AllowsForSpacesAndTabsInFirstLine) {
  std::stringstream ss("  \t2 \t   \n\n"
                       "H      0.0  0.0  -2.0\n"
                       "V      3.0  0.4   0.0\n");

  ASSERT_NO_THROW(XyzStreamHandler::read(ss));
}

TEST_F(XyzStreamHandlerTest, ThrowsForInvalidLineInput) {
  std::stringstream ss1("2\n\n"
                        "H      0.0  0k0  -2.0\n"
                        "V      3.0  0.4   0.0\n");
  std::stringstream ss2("2\n\n"
                        "H      0.0  0k0  -2.0\n"
                        "V      3.0  0.4      \n");
  std::stringstream ss3("2\n\n"
                        "5      0.0  0k0  -2.0\n"
                        "V      3.0  0.4      \n");

  ASSERT_THROW(XyzStreamHandler::read(ss1), FormattedStreamHandler::FormatMismatchException);
  ASSERT_THROW(XyzStreamHandler::read(ss2), FormattedStreamHandler::FormatMismatchException);
  ASSERT_THROW(XyzStreamHandler::read(ss3), FormattedStreamHandler::FormatMismatchException);
}

TEST_F(XyzStreamHandlerTest, AllowsForFurtherColumns) {
  std::stringstream ss1("2\n\n"
                        "H      0.0  0.0  -2.0 and some text later on\n"
                        "V      3.0  0.4   0.0\n");
  std::stringstream ss2("2\n\n"
                        "H      0.0  0.0  -2.0    0.4\n"
                        "V      3.0  0.4   0.0    0.4\n");

  ASSERT_NO_THROW(XyzStreamHandler::read(ss1));
  ASSERT_NO_THROW(XyzStreamHandler::read(ss2));
}

TEST_F(XyzStreamHandlerTest, ThrowsForInvalidAtomTypes) {
  std::stringstream ss1("2\n\n"
                        "Zu     0.0  0.0  -2.0\n"
                        "Vul    3.0  0.4   0.0\n");

  ASSERT_THROW(XyzStreamHandler::read(ss1), FormattedStreamHandler::FormatMismatchException);
}

TEST_F(XyzStreamHandlerTest, AllowsLeadingSpaces) {
  // F.i. as in ORCA-generated xyz files.
  std::stringstream ss("2\n\n"
                       " H      0.0  0.0  -2.0\n"
                       "   V      3.0  0.4   0.0\n");

  auto structure = XyzStreamHandler::read(ss);

  auto expectedElements = ElementTypeCollection{ElementType::H, ElementType::V};
  PositionCollection expectedPositions(2, 3);
  expectedPositions << 0, 0, -2 * Constants::bohr_per_angstrom, 3 * Constants::bohr_per_angstrom,
      0.4 * Constants::bohr_per_angstrom, 0;

  ASSERT_THAT(structure.size(), Eq(2));
  ASSERT_THAT(structure.getElements(), Eq(expectedElements));
  ASSERT_THAT(structure.getPositions(), Eq(expectedPositions));
}

TEST_F(XyzStreamHandlerTest, AllowsForAllUppercaseSymbols) {
  // F.i. as generated in some packages
  std::stringstream ss("2\n\n"
                       "HE      0.0  0.0  -2.0\n"
                       "CA      3.0  0.4   0.0\n");

  auto structure = XyzStreamHandler::read(ss);

  auto expectedElements = ElementTypeCollection{ElementType::He, ElementType::Ca};

  ASSERT_THAT(structure.size(), Eq(2));
  ASSERT_THAT(structure.getElements(), Eq(expectedElements));
}

TEST_F(XyzStreamHandlerTest, AllowsForAllLowercaseSymbols) {
  // F.i. as generated in some packages
  std::stringstream ss("2\n\n"
                       "he      0.0  0.0  -2.0\n"
                       "ca      3.0  0.4   0.0\n");

  auto structure = XyzStreamHandler::read(ss);

  auto expectedElements = ElementTypeCollection{ElementType::He, ElementType::Ca};

  ASSERT_THAT(structure.size(), Eq(2));
  ASSERT_THAT(structure.getElements(), Eq(expectedElements));
}

TEST_F(XyzStreamHandlerTest, CorrectImportInLocaleWithCommas) {
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

    auto structure = XyzStreamHandler::read(ss);

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

TEST_F(XyzStreamHandlerTest, SelfConsistent) {
  std::stringstream in(dimethylbutaneXYZ), out;
  AtomCollection atoms = XyzStreamHandler::read(in);
  XyzStreamHandler::write(out, atoms);
  AtomCollection secondTime = XyzStreamHandler::read(out);

  // Compare AtomCollections
  ASSERT_THAT(atoms, Eq(secondTime));
}

TEST_F(XyzStreamHandlerTest, CommentIsCorrect) {
  std::stringstream in(dimethylbutaneXYZ), out;
  AtomCollection atoms = XyzStreamHandler::read(in);
  XyzStreamHandler::write(out, atoms, "comment");
  std::string line;
  std::istringstream f(out.str());
  /* Ignore first line */
  std::getline(f, line);
  std::getline(f, line);
  ASSERT_THAT(line, Eq("comment"));
}

TEST_F(XyzStreamHandlerTest, CanReadNuclearElectronic) {
  std::stringstream ss("2\n\n"
                       "H      0.0  0.0  -2.0 q\n"
                       "V      3.0  0.4   0.0 Q\n");

  auto structure = XyzStreamHandler::readNuclearElectronic(ss);

  auto expectedElements = ElementTypeCollection{ElementType::H, ElementType::V};
  PositionCollection expectedPositions(2, 3);
  expectedPositions << 0, 0, -2 * Constants::bohr_per_angstrom, 3 * Constants::bohr_per_angstrom,
      0.4 * Constants::bohr_per_angstrom, 0;
  ASSERT_THAT(structure.first.size(), Eq(2));
  ASSERT_THAT(structure.first.getElements(), Eq(expectedElements));
  ASSERT_THAT(structure.first.getPositions(), Eq(expectedPositions));
}

TEST_F(XyzStreamHandlerTest, AsignQCorrectlyNuclearElectronic) {
  std::stringstream ss("3\n\n"
                       "H      0.0  0.0  -2.0 q\n"
                       "O      0.0  0.0   3.0 \n"
                       "V      3.0  0.4   0.0 Q\n");

  auto structure = XyzStreamHandler::readNuclearElectronic(ss);

  ASSERT_EQ(structure.second[0], true);
  ASSERT_EQ(structure.second[1], false);
  ASSERT_EQ(structure.second[2], true);
}

TEST_F(XyzStreamHandlerTest, ReadFileWithAbsoluteLargeCoordinates) {
  const boost::filesystem::path pathToResource =
      boost::dll::program_location().parent_path() / "Resources/large_coordinates.xyz";
  std::ifstream input(pathToResource);
  auto atoms = XyzStreamHandler::read(input);
  input.close();
  std::string temporaryFilePath = "temporaryXYZFile.xyz";
  std::ofstream output;
  output.open(temporaryFilePath);
  EXPECT_TRUE(output.is_open());
  XyzStreamHandler::write(output, atoms);
  output.close();
  std::ifstream input2(temporaryFilePath);
  EXPECT_TRUE(input2.is_open());
  EXPECT_NO_THROW(XyzStreamHandler::read(input2));
  std::remove(temporaryFilePath.c_str());
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
