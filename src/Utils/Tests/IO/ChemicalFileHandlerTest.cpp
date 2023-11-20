/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/Constants.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <gmock/gmock.h>
#include <boost/filesystem.hpp>
#include <fstream>

using namespace testing;
namespace Scine {
namespace Utils {
namespace Tests {

// These are provided by MOLStreamIOTest
extern std::string asymCarbon;
extern std::string dimethylbutane;
// This is provided by XYZStreamIOTest
extern std::string dimethylbutaneXYZ;

TEST(ChemicalFileHandler, SelfConsistentMOLAtomCollectionAndBondOrders) {
  const std::string file = "ChemicalFileHandlerTestArtifact.mol";

  for (const auto& MOLInput : {asymCarbon, dimethylbutane}) {
    std::ofstream out(file);
    out << MOLInput;
    out.close();

    auto fullData = ChemicalFileHandler::read(file);
    ChemicalFileHandler::write(file, fullData.first, fullData.second);
    auto rereadData = ChemicalFileHandler::read(file);

    ASSERT_THAT(fullData, Eq(rereadData));
  }

  boost::filesystem::remove(file);
}

TEST(ChemicalFileHandler, SelfConsistentXYZAtomCollection) {
  const std::string file = "ChemicalFileHandlerTestArtifact.xyz";

  for (const auto& XYZInput : {dimethylbutaneXYZ}) {
    std::ofstream out(file);
    out << XYZInput;
    out.close();

    auto allData = ChemicalFileHandler::read(file);
    ChemicalFileHandler::write(file, allData.first);
    auto rereadData = ChemicalFileHandler::read(file);

    ASSERT_THAT(allData.first, Eq(rereadData.first));
  }

  boost::filesystem::remove(file);
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
