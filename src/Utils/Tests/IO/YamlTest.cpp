/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <gmock/gmock.h>

using namespace testing;

#include <Utils/IO/Yaml.h>
#include <Utils/UniversalSettings/ParametrizedOptionValue.h>
#include <yaml-cpp/yaml.h>

namespace Scine {
namespace Utils {
namespace Tests {

TEST(Yaml, DeserializeValueCollection) {
  /* Make a complicated value collection */
  UniversalSettings::ValueCollection fullCollection;
  fullCollection.addBool("bool", false);
  fullCollection.addInt("int", 4);
  fullCollection.addDouble("float", 0.5);
  fullCollection.addString("str", "Green eggs and ham");
  fullCollection.addIntList("intlist", std::vector<int>{{1, 2, 3}});
  fullCollection.addDoubleList("doublelist", std::vector<double>{{4.0, 5.0, 6.0}});
  fullCollection.addStringList("strlist", std::vector<std::string>{{"hello", "there", "world"}});

  std::vector<UniversalSettings::ValueCollection> collectionList;
  UniversalSettings::ValueCollection firstFruitBowl;
  firstFruitBowl.addInt("apples", 4);
  firstFruitBowl.addInt("bananas", 3);
  collectionList.push_back(firstFruitBowl);

  UniversalSettings::ValueCollection secondFruitBowl;
  secondFruitBowl.addInt("cherries", 10);
  secondFruitBowl.addInt("dates", 0);
  collectionList.push_back(secondFruitBowl);

  fullCollection.addCollectionList("colllist", collectionList);

  UniversalSettings::ParametrizedOptionValue option{"second_fruit_bowl", secondFruitBowl};
  fullCollection.addOptionWithSettings("option", option);

  UniversalSettings::ValueCollection kitchen;
  kitchen.addString("garbage", "full");
  kitchen.addBool("fridge", true);
  fullCollection.addCollection("kitchen", kitchen);

  /* Central symmetry test */
  ASSERT_EQ(fullCollection, deserializeValueCollection(YAML::Load(yamlSerialize(fullCollection))));
}

} /* namespace Tests */
} /* namespace Utils */
} /* namespace Scine */
