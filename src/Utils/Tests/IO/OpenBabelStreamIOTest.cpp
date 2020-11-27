/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/Constants.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/OpenBabelStreamHandler.h>
#include <Utils/Technical/ScopedLocale.h>
#include <gmock/gmock.h>
#include <boost/optional.hpp>
#include <sstream>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

// Provided by MOLStreamIOTest
extern std::string asymCarbon;
extern std::string dimethylbutane;

TEST(OpenBabelStreamIO, FormatListIsOrdered) {
  if (OpenBabelStreamHandler::checkForBinary()) {
    OpenBabelStreamHandler handler;

    const auto& formatList = handler.getSupportedFormats();
    ASSERT_TRUE(std::is_sorted(std::begin(formatList), std::end(formatList)));

    ASSERT_TRUE(handler.formatSupported("mol", OpenBabelStreamHandler::SupportType::ReadWrite));
    ASSERT_TRUE(handler.formatSupported("can", OpenBabelStreamHandler::SupportType::ReadWrite));
    ASSERT_TRUE(handler.formatSupported("smi", OpenBabelStreamHandler::SupportType::ReadWrite));
    ASSERT_TRUE(handler.formatSupported("smiles", OpenBabelStreamHandler::SupportType::ReadWrite));
    ASSERT_TRUE(handler.formatSupported("inchi", OpenBabelStreamHandler::SupportType::ReadWrite));
  }
}

TEST(OpenBabelStreamIO, ReadWriteSMILES) {
  std::vector<std::string> exampleSMILES{{
      "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",             // Caffeine
      "CC(=O)C1CCC2C1(CCC3C2CCC4=CC(=O)CCC34C)C", // Progesterone, canonical
      "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C"       // Testosterone, canonical
  }};

  if (OpenBabelStreamHandler::checkForBinary()) {
    OpenBabelStreamHandler handler;
    for (const std::string& SMILESInput : exampleSMILES) {
      std::stringstream in(SMILESInput), out;
      std::pair<AtomCollection, BondOrderCollection> allData;
      ASSERT_NO_THROW(allData = handler.read(in, "can"));
      ASSERT_NO_THROW(handler.write(out, "can", allData.first, allData.second));
    }
  }
}

TEST(OpenBabelStreamIO, ReadWriteInChI) {
  std::vector<std::string> exampleInChI = {
      std::string(R"(InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3)"), // Caffeine
      std::string(R"(InChI=1S/C21H30O2/c1-13(22)17-6-7-18-16-5-4-14-12-15(23)8-10-20(14,2)19(16)9-11-21(17,18)3/h12,16-19H,4-11H2,1-3H3/t16-,17+,18-,19-,20-,21+/m0/s1)") // Progesterone
      /* For some reason, trying testosterone InChI will make obabel fall flat on
       * its face with a wacky mol output. See for yourself:
       *
       * 1. $ obabel -iinchi -omol --gen3D -h
       * 2. copy in testosterone's InChI string
       * 3. Press ctrl-d
       */
  };

  if (OpenBabelStreamHandler::checkForBinary()) {
    OpenBabelStreamHandler handler;
    for (const std::string& InChIInput : exampleInChI) {
      std::stringstream in(InChIInput), out;
      std::pair<AtomCollection, BondOrderCollection> allData;
      ASSERT_NO_THROW(allData = handler.read(in, "inchi"));
      ASSERT_NO_THROW(handler.write(out, "inchi", allData.first, allData.second));
    }
  }
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
