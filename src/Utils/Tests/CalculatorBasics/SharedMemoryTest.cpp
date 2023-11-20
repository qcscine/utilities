/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Core/Log.h>
#include <Utils/CalculatorBasics.h>
#include <gmock/gmock.h>
#include <unistd.h>
#include <chrono>
#include <thread>

using namespace testing;

namespace Scine {
namespace Utils {
namespace Tests {

class ASharedMemoryTest : public Test {
 public:
  std::shared_ptr<SharedMemory::MemoryManager> manager;

 private:
  void SetUp() final {
    auto log = Core::Log::silent();
    manager = std::make_shared<SharedMemory::MemoryManager>(log);
  }
};

TEST_F(ASharedMemoryTest, CanCommunicate) {
  // prepare long random string to ensure overflow safety
  static const char alphanum[] = "0123456789"
                                 "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                 "abcdefghijklmnopqrstuvwxyz";
  int len = 1000;
  std::string longString;
  longString.reserve(len);
  for (int i = 0; i < len; ++i) {
    longString += alphanum[rand() % (sizeof(alphanum) - 1)];
  }
  manager->setStatus(SharedMemory::status::notStarted);
  int pid = fork();
  if (pid == -1) {
    throw std::runtime_error("Failed to fork");
  }
  else if (pid == 0) {
    manager->setStatus(SharedMemory::status::running);

    while (manager->getStatus() == SharedMemory::status::running) {
      std::this_thread::sleep_for(std::chrono::microseconds(10));
    }

    manager->setEnergy(42.0);
    ASSERT_TRUE(std::fabs(manager->getEnergy() - 42.0) < 1e-6);
    manager->setStatus(SharedMemory::status::running);

    while (manager->getStatus() == SharedMemory::status::running) {
      std::this_thread::sleep_for(std::chrono::microseconds(10));
    }

    manager->setError("some message");
    ASSERT_TRUE(manager->getError() == "some message");
    manager->setStatus(SharedMemory::status::running);

    while (manager->getStatus() == SharedMemory::status::running) {
      std::this_thread::sleep_for(std::chrono::microseconds(10));
    }

    manager->setError(longString);
    ASSERT_FALSE(manager->getError().empty());
    ASSERT_TRUE(manager->getError() != "some message");
    manager->setStatus(SharedMemory::status::finished);
    exit(0);
  }
  else {
    while (manager->getStatus() == SharedMemory::status::notStarted) {
      std::this_thread::sleep_for(std::chrono::microseconds(10));
    }

    ASSERT_TRUE(manager->getStatus() == SharedMemory::status::running);
    manager->setStatus(SharedMemory::status::notStarted);

    while (manager->getStatus() == SharedMemory::status::notStarted) {
      std::this_thread::sleep_for(std::chrono::microseconds(10));
    }

    ASSERT_TRUE(std::fabs(manager->getEnergy() - 42.0) < 1e-6);
    manager->setStatus(SharedMemory::status::notStarted);

    while (manager->getStatus() == SharedMemory::status::notStarted) {
      std::this_thread::sleep_for(std::chrono::microseconds(10));
    }

    ASSERT_TRUE(manager->getError() == "some message");
    manager->setError("");
    manager->setStatus(SharedMemory::status::notStarted);

    while (manager->getStatus() == SharedMemory::status::notStarted) {
      std::this_thread::sleep_for(std::chrono::microseconds(10));
    }
    ASSERT_FALSE(manager->getError().empty());
    ASSERT_TRUE(manager->getStatus() == SharedMemory::status::finished);
  }
  ASSERT_FALSE(manager->wasCleanedUp());
  manager->cleanUp();
  ASSERT_TRUE(manager->wasCleanedUp());
}

} // namespace Tests
} // namespace Utils
} // namespace Scine
