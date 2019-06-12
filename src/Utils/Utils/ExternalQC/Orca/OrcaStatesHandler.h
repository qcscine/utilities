/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_ORCASTATESHANDLER_H
#define UTILS_ORCASTATESHANDLER_H

#include <Utils/CalculatorBasics/StatesHandler.h>

namespace Scine {
namespace Utils {
namespace ExternalQC {
class OrcaCalculator;

class StateNotCompatibleException : public std::exception {
 public:
  const char* what() const noexcept final {
    return "State is not compatible with the ORCA calculator.";
  }
};

class StateSavingException : public std::exception {
 public:
  const char* what() const noexcept final {
    return "ORCA state could not be saved.";
  }
};

class OrcaStatesHandler : public Utils::StatesHandler {
 public:
  /**
   * @brief Constructor, it takes the orca calculator as argument in order to get some information from it.
   * @param orcaCalculator The embedding method wrapper.
   */
  explicit OrcaStatesHandler(OrcaCalculator& calculator);
  ~OrcaStatesHandler() final;
  /**
   * @brief The implementation of the store function. It stores a state with the specified size.
   * @param size The desired Utils::StateSize of the state to store.
   */
  void store(Utils::StateSize size) final;
  /**
   * @brief Loads a state in the ORCA calculator by copying an .gbw file.
   * @param state The state to load. If it is not compatible with the ORCA calculator,
   *              a StateNotCompatibleException is thrown.
   */
  void load(std::shared_ptr<Utils::State> state) final;
  /**
   * @brief Returns the current state.
   * @param The desired Utils::StateSize of the state to store.
   *            This method does not need to store the current state in the stateHandler,
   *            it is a const method and therefore suited to be used in a copy-constructor.
   */
  std::shared_ptr<Utils::State> getCurrentState(Utils::StateSize size) const final;

 private:
  // Copying .gbw files
  void copyBackupFile(const std::string& from, const std::string& to);
  OrcaCalculator& calculator_;
};

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine

#endif // UTILS_ORCASTATESHANDLER_H
