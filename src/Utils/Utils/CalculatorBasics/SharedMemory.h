/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher
 * Group.\n See LICENSE.txt for details.
 */
#ifndef UTILS_SHAREDMEMORY_H_
#define UTILS_SHAREDMEMORY_H_

#include <Core/Log.h>
#include <sys/ipc.h>   /* general SysV IPC structures          */
#include <sys/sem.h>   /* shared memory functions and structs. */
#include <sys/shm.h>   /* semaphore functions and structs.     */
#include <sys/types.h> /* various type definitions.            */
#include <cstdio>
#include <string>

namespace Scine {

namespace Utils {

namespace SharedMemory {

/**
 * @enum status SharedMemory.h
 * @brief An enum holding different possible status to be held in shared memory
 */
enum status { failed = -1, notStarted = 0, running = 1, finished = 2, aborted = 3 };

#define MAX_ERROR 496

/**
 * @class Memory SharedMemory.h
 * @brief A struct containing everything that has to be communicated between forks
 */
struct Memory {
  /// status_code in shared mem to communicate
  int status = 0;
  // energy value
  double energy = 0.0;
  /// error message
  char error[MAX_ERROR] = "";
};

/**
 * @class MemoryManager SharedMemory.h
 * @brief A class handling access and holding the shared memory.
 *
 * inspired by approach by
 * http://www.cs.kent.edu/~ruttan/sysprog/lectures/shmem/shared-mem-with-semaphore.c
 * (accessed August 2022)
 * but packed into a class and hide all the direct memory access
 * memory is allocated and deallocated in constructor and cleanUp method
 * memory interaction + semaphore lock is only possible via getter/setter
 */
class MemoryManager {
 public:
  explicit MemoryManager(Core::Log& log) {
    _log = std::make_shared<Core::Log>(log);
    // set up shared memory
    // in case of debugging, `perror()` can be called before throws
    // code has been minimally altered in comparison to the given reference
    // above

    /* create a semaphore set with ID 250, with one semaphore   */
    /* in it, with access only to the owner.                    */
    semSetId = semget(250, 1, IPC_CREAT | 0600);
    if (semSetId == -1) {
      perror("");
      throw std::runtime_error("Failed to create semaphore set for memory management");
    }
    /* initialize the first (and single) semaphore in our set to '1'. */
    semVal.val = 1;
    int rc = semctl(semSetId, 0, SETVAL, semVal);
    if (rc == -1) {
      perror("");
      throw std::runtime_error("Failed to initialize a semaphore for memory management");
    }
    /* allocate a shared memory segment */
    shmId = shmget(100, sizeof(struct Memory), IPC_CREAT | 0600);
    if (shmId == -1) {
      perror("");
      throw std::runtime_error("Failed to allocate shared memory");
    }
    /* attach the shared memory segment to our process's address space. */
    shmAddr = (char*)shmat(shmId, NULL, 0);
    if (!shmAddr) { /* operation failed. */
      perror("");
      throw std::runtime_error("Failed to attach the shared memory segment to "
                               "our process's address space.");
    }
    /* set the memory pointer to the shared memory segment. */
    memory = (struct Memory*)shmAddr;
  }
  ~MemoryManager() {
    if (!_wasCleanedUp) {
      try {
        cleanUp();
      }
      catch (...) {
        // here we catch everything, because destructors must not throw exceptions
      }
    }
  }

  /**
   * @brief Needs to be called after parallel session to clean up shared memory
   */
  inline void cleanUp() {
    _wasCleanedUp = true;
    /* de-allocate the shared memory segment. */
    if (shmctl(shmId, IPC_RMID, NULL) == -1) {
      perror("");
      throw std::runtime_error("Failed to deallocate the shared memory segment");
    }
    /* detach the shared memory segment from our process's address space. */
    if (shmdt(shmAddr) == -1) {
      perror("");
      throw std::runtime_error("Failed to detach the shared memory segment");
    }
  };

  /**
   * @brief Get the status of the calculation from the shared memory with
   * included semaphore locking
   */
  inline status getStatus() const {
    if (_wasCleanedUp) {
      throw std::runtime_error("Shared memory has already been cleaned up");
    }
    semLock(semSetId);
    auto status = static_cast<SharedMemory::status>(memory->status);
    semUnlock(semSetId);
    return status;
  }
  /**
   * @brief set the status of the calculation in the shared memory with included
   * semaphore locking
   * @param s The status
   */
  inline void setStatus(status s) {
    if (_wasCleanedUp) {
      throw std::runtime_error("Shared memory has already been cleaned up");
    }
    semLock(semSetId);
    memory->status = static_cast<int>(s);
    semUnlock(semSetId);
  }
  /**
   * @brief Get the energy of the calculation from the shared memory with
   * included semaphore locking
   */
  inline double getEnergy() const {
    if (_wasCleanedUp) {
      throw std::runtime_error("Shared memory has already been cleaned up");
    }
    semLock(semSetId);
    double e = memory->energy;
    semUnlock(semSetId);
    return e;
  }
  /**
   * @brief set the energy of the calculation in the shared memory with included
   * semaphore locking
   * @param e The energy value
   */
  inline void setEnergy(double e) {
    if (_wasCleanedUp) {
      throw std::runtime_error("Shared memory has already been cleaned up");
    }
    semLock(semSetId);
    memory->energy = e;
    semUnlock(semSetId);
  }
  /**
   * @brief Get the error message from the shared memory with
   * included semaphore locking
   */
  inline std::string getError() const {
    if (_wasCleanedUp) {
      throw std::runtime_error("Shared memory has already been cleaned up");
    }
    semLock(semSetId);
    auto error = std::string(memory->error);
    semUnlock(semSetId);
    return error;
  }
  /**
   * @brief set the error message in the shared memory with included
   * semaphore locking
   * @param s The error message
   */
  inline void setError(std::string s) {
    if (_wasCleanedUp) {
      throw std::runtime_error("Shared memory has already been cleaned up");
    }
    semLock(semSetId);
    if (s.size() > MAX_ERROR) {
      _log->error << "Received too large error message for sharing between processes:\n" << s << Core::Log::endl;
      s = s.substr(0, MAX_ERROR);
    }
    strcpy(memory->error, s.c_str());
    semUnlock(semSetId);
  }
  /**
   * @brief return if shared memory has been cleaned up already or not
   */
  inline bool wasCleanedUp() const {
    return _wasCleanedUp;
  }

 private:
  /**
   * @brief locks the semaphore, for exclusive access to a resource
   * @param semSetId The semaphore set ID
   */
  inline void semLock(int semSetId) const {
    /* structure for semaphore operations.   */
    struct sembuf semOp;
    /* wait on the semaphore, unless it's value is non-negative. */
    semOp.sem_num = 0;
    semOp.sem_op = -1;
    semOp.sem_flg = 0;
    semop(semSetId, &semOp, 1);
  }

  /**
   * @brief unlocks the semaphore
   * @param semSetId The semaphore set ID
   */
  inline void semUnlock(int semSetId) const {
    /* structure for semaphore operations.   */
    struct sembuf semOp;
    /* signal the semaphore - increase its value by one. */
    semOp.sem_num = 0;
    semOp.sem_op = 1;
    semOp.sem_flg = 0;
    semop(semSetId, &semOp, 1);
  }
  /// ID of the semaphore set
  int semSetId;
  /// semaphore value, for semctl()
  union semun {
    int val;
    struct semid_ds* buf;
    ushort* array;
  } semVal;
  /// ID of the shared memory segment
  int shmId;
  /// address of shared memory segment
  char* shmAddr;
  /// the object in shared memory containing information
  struct Memory* memory;
  std::shared_ptr<Core::Log> _log;
  bool _wasCleanedUp = false;
};

} // namespace SharedMemory
} // namespace Utils
} // namespace Scine
#endif // UTILS_SHAREDMEMORY_H_
