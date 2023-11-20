/**
 * @file CloneInterface.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Technical/CloneInterface.h>
#include <gmock/gmock.h>
#include <memory>

using namespace testing;

namespace Scine {
namespace Utils {
namespace Tests {

struct TestObject {
  TestObject(double d) : testDouble(d) {
  }
  double testDouble;
};

class InterfaceTest {
 public:
  virtual ~InterfaceTest() = default;
  std::shared_ptr<InterfaceTest> clone() const {
    return this->cloneImpl();
  }

  virtual double getNumber() const noexcept = 0;

 private:
  virtual std::shared_ptr<InterfaceTest> cloneImpl() const = 0;
};

class AbstractTest : public CloneInterface<Abstract<AbstractTest>, InterfaceTest> {
 public:
  virtual void addOne() noexcept = 0;
};

class ConcreteNoAbstract : public CloneInterface<ConcreteNoAbstract, InterfaceTest> {
 public:
  ConcreteNoAbstract(double d) {
    testPtr = std::make_unique<TestObject>(d);
  }
  ConcreteNoAbstract(const ConcreteNoAbstract& rhs) : ConcreteNoAbstract(rhs.getNumber()) {
  }

  double getNumber() const noexcept final {
    return testPtr->testDouble;
  }

  std::unique_ptr<TestObject> testPtr;
};

class ConcreteWithAbstract : public CloneInterface<ConcreteWithAbstract, AbstractTest, InterfaceTest> {
 public:
  ConcreteWithAbstract(double d) {
    testPtr = std::make_unique<TestObject>(d);
  }
  ConcreteWithAbstract(const ConcreteWithAbstract& rhs) : ConcreteWithAbstract(rhs.getNumber()) {
  }

  double getNumber() const noexcept final {
    return testPtr->testDouble;
  }

  void addOne() noexcept final {
    ++(testPtr->testDouble);
  }

  void addTwo() noexcept {
    testPtr->testDouble += 2;
  }

  std::unique_ptr<TestObject> testPtr;
};

class ACloneInterfaceTest : public Test {
 public:
  std::unique_ptr<InterfaceTest> interface;
  std::unique_ptr<AbstractTest> abstract;
  std::unique_ptr<ConcreteWithAbstract> concrete;
};

TEST_F(ACloneInterfaceTest, CanCloneConcreteNoAbstractClass) {
  interface = std::make_unique<ConcreteNoAbstract>(2.0);
  auto clonedInterface = interface->clone();
  interface = std::make_unique<ConcreteNoAbstract>(3.0);
  auto clonedInterface3 = interface->clone();

  ASSERT_EQ(clonedInterface->getNumber(), 2.0);
  ASSERT_EQ(clonedInterface3->getNumber(), 3.0);
}

TEST_F(ACloneInterfaceTest, CanCloneConcreteWithAbstractClass) {
  interface = std::make_unique<ConcreteWithAbstract>(2.0);
  auto clonedInterface = interface->clone();

  interface = std::make_unique<ConcreteWithAbstract>(3.0);
  auto clonedInterface3 = interface->clone();

  abstract = std::make_unique<ConcreteWithAbstract>(2.0);
  abstract->addOne();
  auto clonedFromAbstract = abstract->clone();

  ASSERT_EQ(clonedInterface->getNumber(), 2.0);
  ASSERT_EQ(clonedInterface3->getNumber(), 3.0);
  ASSERT_EQ(clonedFromAbstract->getNumber(), 3.0);
}

TEST_F(ACloneInterfaceTest, CloneMethodHidesBaseClassMethod) {
  concrete = std::make_unique<ConcreteWithAbstract>(2.0);
  auto clonedConcrete = concrete->clone();

  clonedConcrete->addOne();
  clonedConcrete->addTwo();

  ASSERT_EQ(clonedConcrete->getNumber(), 5.0);
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
