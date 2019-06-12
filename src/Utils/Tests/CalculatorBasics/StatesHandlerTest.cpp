/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/CalculatorBasics/State.h>
#include <Utils/CalculatorBasics/StatesHandler.h>
#include <gmock/gmock.h>

namespace Scine {
namespace Utils {
namespace Tests {

using namespace testing;

class AStatesHandlerTest : public Test {
 public:
  // Mock State class for testing the StateHandler.
  class TestState : public State {
   public:
    explicit TestState(StateSize size) : State(size) {
    }
    const Eigen::MatrixXd& getMatrixState(const std::string& matrixState) const final {
      static Eigen::MatrixXd f;
      return f;
    };
    const std::string& getStringState(const std::string& stringState) const final {
      static std::string f;
      return f;
    };
    int getIntState(const std::string& integerState) const final {
      return 0;
    };
    double getDoubleState(const std::string& doubleState) const final {
      return 0.0;
    };
    void initialize() final{};
  };

  // Mock StatesHandler class for testing the StatesHandler.
  class TestStatesHandler : public StatesHandler {
   public:
    ~TestStatesHandler() final = default;
    void store(StateSize size) final {
      TestState state(size);
      states_.emplace_back(std::make_unique<TestState>(state));
    };
    std::shared_ptr<State> getCurrentState(StateSize /*size*/) const final {
      return std::make_shared<TestState>(TestState(StateSize::minimal));
    }
  };

  TestStatesHandler statesHandler_;
};

TEST_F(AStatesHandlerTest, CanStoreTestState) {
  statesHandler_.store(StateSize::minimal);
  statesHandler_.store(StateSize::regular);
  statesHandler_.store(StateSize::extensive);
}

TEST_F(AStatesHandlerTest, CanRecoverTestState) {
  statesHandler_.store(StateSize::minimal);
  statesHandler_.store(StateSize::regular);
  statesHandler_.store(StateSize::extensive);

  auto newestState = statesHandler_.popNewestState();
  auto middleState = statesHandler_.popNewestState();
  auto oldestState = statesHandler_.popNewestState();

  ASSERT_EQ(newestState->getStateSize(), StateSize::extensive);
  ASSERT_EQ(middleState->getStateSize(), StateSize::regular);
  ASSERT_EQ(oldestState->getStateSize(), StateSize::minimal);
}

TEST_F(AStatesHandlerTest, CanRecoverTestStateFromOldest) {
  statesHandler_.store(StateSize::minimal);
  statesHandler_.store(StateSize::regular);
  statesHandler_.store(StateSize::extensive);

  auto oldestState = statesHandler_.popOldestState();
  auto middleState = statesHandler_.popOldestState();
  auto newestState = statesHandler_.popOldestState();

  ASSERT_EQ(newestState->getStateSize(), StateSize::extensive);
  ASSERT_EQ(middleState->getStateSize(), StateSize::regular);
  ASSERT_EQ(oldestState->getStateSize(), StateSize::minimal);
}

TEST_F(AStatesHandlerTest, ThrowsIfNoStatesLeftButOneRequired) {
  statesHandler_.store(StateSize::minimal);
  statesHandler_.popNewestState();

  ASSERT_THROW(statesHandler_.popNewestState(), EmptyStatesHandlerContainer);
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
