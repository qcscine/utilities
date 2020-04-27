/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/Constants.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/MolStreamHandler.h>
#include <Utils/Technical/ScopedLocale.h>
#include <gmock/gmock.h>
#include <boost/optional.hpp>
#include <sstream>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {
/**
 * @class MolStreamHandlerTest
 * @brief Comprises tests for the class Scine::Utils::AtomCollection.
 * @test
 */
class MolStreamHandlerTest : public Test {
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

std::string asymCarbon = R"delim(57424940
-OEChem-11011611493D

  5  4  0     1  0  0  0  0  0999 V2000
    1.8992   -0.5423    0.3511 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2940   -1.3226    0.4558 Br  0  0  0  0  0  0  0  0  0  0  0  0
   -0.4101    1.5638    0.7393 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -0.1538    0.2476   -1.4431 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0413    0.0536   -0.1030 C   0  0  1  0  0  0  0  0  0  0  0  0
  1  5  1  0  0  0  0
  2  5  1  0  0  0  0
  3  5  1  0  0  0  0
  4  5  1  0  0  0  0
M  ISO  1  1  2
M  END)delim";

std::string dimethylbutane = "C6H14 Butane, 2,2-dimethyl- 75832\n"
                             "##CCCBDB10101603:55\n"
                             "Geometry Optimized at B3LYP/TZVP\n"
                             " 20 19  0  0  0  0  0  0  0  0999 V2000\n"
                             "    1.9046   -1.0745    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                             "    1.1702    0.2696    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                             "   -0.3765    0.2200    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                             "   -0.8983    1.6682    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                             "   -0.8983   -0.4979    1.2565 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                             "   -0.8983   -0.4979   -1.2565 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                             "   -1.9910    1.6914    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                             "    2.9849   -0.9134    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                             "   -0.5247   -0.0194    2.1656 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                             "   -0.5247   -0.0194   -2.1656 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                             "   -0.5966   -1.5471    1.2808 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                             "   -0.5966   -1.5471   -1.2808 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                             "   -1.9906   -0.4711    1.2894 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                             "   -1.9906   -0.4711   -1.2894 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                             "    1.6651   -1.6718    0.8821 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                             "    1.6651   -1.6718   -0.8821 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                             "    1.4911    0.8456    0.8748 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                             "    1.4911    0.8456   -0.8748 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                             "   -0.5520    2.2123    0.8829 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                             "   -0.5520    2.2123   -0.8829 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
                             "  1  2  1  0  0  0\n"
                             "  1  8  1  0  0  0\n"
                             "  1 15  1  0  0  0\n"
                             "  1 16  1  0  0  0\n"
                             "  2  3  1  0  0  0\n"
                             "  2 17  1  0  0  0\n"
                             "  2 18  1  0  0  0\n"
                             "  3  4  1  0  0  0\n"
                             "  3  5  1  0  0  0\n"
                             "  3  6  1  0  0  0\n"
                             "  4  7  1  0  0  0\n"
                             "  4 19  1  0  0  0\n"
                             "  4 20  1  0  0  0\n"
                             "  5  9  1  0  0  0\n"
                             "  5 11  1  0  0  0\n"
                             "  5 13  1  0  0  0\n"
                             "  6 10  1  0  0  0\n"
                             "  6 12  1  0  0  0\n"
                             "  6 14  1  0  0  0\n"
                             "M  END";

TEST_F(MolStreamHandlerTest, SelfConsistent) {
  for (const auto& MOLInput : {asymCarbon, dimethylbutane}) {
    std::stringstream in(MOLInput), out;
    auto allData = MolStreamHandler::read(in);
    MolStreamHandler::write(out, allData.first, allData.second, "V2000");
    auto secondTime = MolStreamHandler::read(out);

    // Compare AtomCollections
    ASSERT_THAT(allData.first, Eq(secondTime.first));
    // Compare BondOrderCollections
    ASSERT_THAT(allData.second, Eq(secondTime.second));
  }
}

TEST_F(MolStreamHandlerTest, CorrectImportInLocaleWithCommas) {
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
    std::stringstream ss(dimethylbutane);

    auto data = MolStreamHandler::read(ss);

    auto expectedElements = ElementTypeCollection{
        ElementType::C, ElementType::C, ElementType::C, ElementType::C, ElementType::C, ElementType::C, ElementType::H,
        ElementType::H, ElementType::H, ElementType::H, ElementType::H, ElementType::H, ElementType::H, ElementType::H,
        ElementType::H, ElementType::H, ElementType::H, ElementType::H, ElementType::H, ElementType::H};

    ASSERT_THAT(data.first.size(), Eq(expectedElements.size()));
    ASSERT_THAT(data.first.getElements(), Eq(expectedElements));
  }
  catch (std::runtime_error& e) {
    // Do nothing.
  }
}
} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
