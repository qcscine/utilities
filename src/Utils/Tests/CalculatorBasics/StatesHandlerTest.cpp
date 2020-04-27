/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/CalculatorBasics/StatesHandler.h"
#include <Core/BaseClasses/StateHandableObject.h>
#include <Core/Exceptions.h>
#include <gmock/gmock.h>

namespace Scine {
namespace Utils {
namespace Tests {

using namespace testing;

class AStatesHandlerTest : public Test {
 public:
  // Mock State class for testing the StateHandler.
  class TestState : public Core::State {
   public:
    TestState() = default;
    ~TestState() final = default;
    int stateNumber_ = 0;
  };

  class TestStateWithWrongType : public Core::State {
   public:
    TestStateWithWrongType() = default;
    ~TestStateWithWrongType() final = default;
  };

  class TestHandableObject : public Core::StateHandableObject {
   public:
    void loadState(std::shared_ptr<Core::State> state) final {
      auto castState = std::dynamic_pointer_cast<TestState>(state);

      if (!castState)
        throw Core::StateCastingException();
      stateDefiner_ = castState->stateNumber_;
    }

    std::shared_ptr<Core::State> getState() const final {
      TestState state;
      state.stateNumber_ = stateDefiner_;
      return std::make_shared<TestState>(std::move(state));
    }

    void add(int i) {
      stateDefiner_ += i;
    }
    int get() const {
      return stateDefiner_;
    }

   private:
    int stateDefiner_ = 0;
  };

  std::shared_ptr<TestHandableObject> testObject;
  StatesHandler statesHandler_;

  void SetUp() final {
    testObject = std::make_shared<TestHandableObject>();
    statesHandler_ = StatesHandler(testObject);
  }
};

TEST_F(AStatesHandlerTest, CanStoreTestState) {
  statesHandler_.store();
  statesHandler_.store();
  EXPECT_EQ(statesHandler_.size(), 2);
}

TEST_F(AStatesHandlerTest, CanRecoverTestState) {
  statesHandler_.store();
  testObject->add(10);
  statesHandler_.store();
  testObject->add(5);
  EXPECT_EQ(testObject->get(), 15);
  statesHandler_.load(0);
  EXPECT_EQ(testObject->get(), 0);
  statesHandler_.load(1);
  EXPECT_EQ(testObject->get(), 10);
}

TEST_F(AStatesHandlerTest, CanPopOldestTestState) {
  statesHandler_.store();
  testObject->add(10);
  statesHandler_.store();
  testObject->add(5);
  EXPECT_EQ(testObject->get(), 15);
  statesHandler_.load(statesHandler_.popOldestState());
  EXPECT_EQ(testObject->get(), 0);
  EXPECT_EQ(statesHandler_.size(), 1);
}

TEST_F(AStatesHandlerTest, CanPopNewestTestState) {
  statesHandler_.store();
  testObject->add(10);
  statesHandler_.store();
  testObject->add(5);
  EXPECT_EQ(testObject->get(), 15);
  statesHandler_.load(statesHandler_.popNewestState());
  EXPECT_EQ(testObject->get(), 10);
  EXPECT_EQ(statesHandler_.size(), 1);
}
TEST_F(AStatesHandlerTest, ThrowsIfNoStatesLeftButOneRequired) {
  statesHandler_.store();
  statesHandler_.popNewestState();

  ASSERT_THROW(statesHandler_.popNewestState(), EmptyStatesHandlerContainer);
}

TEST_F(AStatesHandlerTest, HandlesMissingInstanceIfExternallyDeleted) {
  testObject = std::make_shared<TestHandableObject>();
  ASSERT_THROW(statesHandler_.store(), NoStateHandableObjectPresent);
}

TEST_F(AStatesHandlerTest, ThrowsIfWrongStateIsLoaded) {
  ASSERT_THROW(statesHandler_.load(std::make_shared<TestStateWithWrongType>()), Core::StateCastingException);
}
} // namespace Tests
} // namespace Utils
} // namespace Scine
