/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/IO/TurbomoleMinimalBasisfile.h>
#include <gmock/gmock.h>
#include <boost/dll/runtime_symbol_info.hpp>
#include <boost/filesystem.hpp>
#include <fstream>

using namespace testing;

namespace Scine {
namespace Utils {
namespace Tests {

class AtomicGtoIOTest : public Test {
 public:
  std::string filepath = "sample.basis";
  std::string atomicGtoBasisfileSample = R"($basis
*
h   STO-6G-PM6NOCORE
*
    6  s
    3.718317372849182e+01    9.834128556209795e-02
    6.817494200625323e+00    1.484283848368730e-01
    1.907289198957793e+00    1.949490512064076e-01
    6.552053163672594e-01    1.923335158882793e-01
    2.544353976087330e-01    1.063406753678467e-01
    1.047905475584809e-01    1.710841115577699e-02
*
i   STO-6G-PM6NOCORE
*
    6  s
    2.853800560769282e+01    2.371960257290400e-02
    1.027654981160127e+01    7.568403684472617e-02
    3.739811981009558e+00   -1.837823897503570e-01
    2.147381809490497e+00   -6.575089050524558e-01
    7.426461452562305e-01    6.201168354362249e-01
    4.479772725495345e-01    1.211347651945825e-01
    6  p
    1.388706457326902e+01    4.445031764907272e-03
    1.285984761665524e+00   -5.701411950542014e-02
    6.185947662045390e-01   -1.079960022734586e-01
    1.986085191388154e-01    1.078398965029914e-01
    1.223430031616701e-01    4.918400151603063e-02
    7.584860834065074e-02    3.416577190190453e-03
    6  d
    3.101543086634682e+00   -8.467148831603270e-02
    1.199346724880346e+00   -9.828743638342213e-02
    3.236493887518321e-01    1.048562083146502e-01
    1.924401591892591e-01    8.125300718126686e-02
    1.192442302993152e-01    1.831114551589139e-02
    7.413120450497401e-02    7.875753473562038e-04
*
$end)";

  void SetUp() override {
    std::ofstream file(filepath);
    file << atomicGtoBasisfileSample;
  }

  void TearDown() override {
    boost::filesystem::remove(filepath);
  }
};

TEST_F(AtomicGtoIOTest, ReadAtomicGtoBasisfile) {
  auto map = readTurbomoleBasisfile(filepath);

  ASSERT_EQ(map.size(), 2);
  ASSERT_TRUE(map.count(1) == 1);
  ASSERT_TRUE(map.at(1).s);
  ASSERT_FALSE(map.at(1).p);
  ASSERT_FALSE(map.at(1).d);
  ASSERT_TRUE(map.count(53) == 1);
  ASSERT_TRUE(map.at(53).s);
  ASSERT_TRUE(map.at(53).p);
  ASSERT_TRUE(map.at(53).d);
}

TEST(STONG, ReadLargeBasisfile) {
  auto pathToResources = boost::dll::program_location().parent_path();
  pathToResources /= "Resources";
  auto map = readTurbomoleBasisfile((pathToResources / "STO-6G.basis").string());
}

TEST_F(AtomicGtoIOTest, AtomicGtoGetNWChemFormat) {
  auto map = readTurbomoleBasisfile(filepath);

  auto nwChemFormatH = map.at(1).getNWChemFormat();
  auto nwChemFormatI = map.at(53).getNWChemFormat();

  ASSERT_EQ(nwChemFormatH.size(), 1);
  ASSERT_EQ(nwChemFormatH[0].size(), 7);

  ASSERT_EQ(nwChemFormatI.size(), 3);
  ASSERT_EQ(nwChemFormatI[0].size(), 7);
  ASSERT_EQ(nwChemFormatI[1].size(), 7);
  ASSERT_EQ(nwChemFormatI[2].size(), 7);

  ASSERT_EQ(boost::get<int>(nwChemFormatH[0][0]), 0);
  auto hPair = boost::get<std::pair<double, double>>(nwChemFormatH[0][2]);
  ASSERT_EQ(hPair.first, 6.817494200625323e+00);
  ASSERT_EQ(hPair.second, 1.484283848368730e-01);

  ASSERT_EQ(boost::get<int>(nwChemFormatI[1][0]), 1);
  ASSERT_EQ(boost::get<int>(nwChemFormatI[2][0]), 2);
  auto iPair = boost::get<std::pair<double, double>>(nwChemFormatI[2][1]);
  ASSERT_EQ(iPair.first, 3.101543086634682e+00);
  ASSERT_EQ(iPair.second, -8.467148831603270e-02);
}

TEST_F(AtomicGtoIOTest, AtomicGtoGetGtfs) {
  auto map = readTurbomoleBasisfile(filepath);

  auto gtfsH = map.at(1).getGtfs();
  auto gtfsI = map.at(53).getGtfs();

  ASSERT_EQ(gtfsH.size(), 1);
  ASSERT_TRUE(gtfsH.find("s") != gtfsH.end());
  ASSERT_FALSE(gtfsH.find("p") != gtfsH.end());
  ASSERT_FALSE(gtfsH.find("d") != gtfsH.end());
  ASSERT_EQ(gtfsH.at("s").size(), 6);

  ASSERT_EQ(gtfsI.size(), 3);
  ASSERT_TRUE(gtfsI.find("s") != gtfsI.end());
  ASSERT_TRUE(gtfsI.find("p") != gtfsI.end());
  ASSERT_TRUE(gtfsI.find("d") != gtfsI.end());
  ASSERT_EQ(gtfsI.at("s").size(), 6);
  ASSERT_EQ(gtfsI.at("p").size(), 6);
  ASSERT_EQ(gtfsI.at("d").size(), 6);

  auto secondSGtfH = gtfsH["s"][1];
  ASSERT_EQ(secondSGtfH.exponent, 6.817494200625323e+00);
  ASSERT_EQ(secondSGtfH.coefficient, 1.484283848368730e-01);

  auto firstDGtfI = gtfsI["d"][0];
  ASSERT_EQ(firstDGtfI.exponent, 3.101543086634682e+00);
  ASSERT_EQ(firstDGtfI.coefficient, -8.467148831603270e-02);
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
