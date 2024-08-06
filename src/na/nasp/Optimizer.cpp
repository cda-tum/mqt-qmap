#include "nasp/Optimizer.hpp"

#include "nasp/Solver.hpp"

#include <algorithm>
#include <cassert>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <functional>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unistd.h>
#include <unordered_map>
#include <utility>
#include <sys/wait.h>
#include <yaml-cpp/yaml.h>

namespace na {

bool childProcess = false;

auto handleAlarm(int) -> void {
  if (childProcess) {
    // kill the process itself
    kill(getpid(), SIGKILL);
  }
}

auto Optimizer::forkChildProcess(
    const std::uint16_t arg, const std::chrono::seconds childTimeout) -> void {
  // Create a pipe
  int pipefd[2];
  if (pipe(pipefd) == -1) {
    killAllChildProcesses();
    throw std::runtime_error("Pipe failed.");
  }
  // Fork the process
  const pid_t pid = fork();
  if (pid < 0) {
    killAllChildProcesses();
    throw std::runtime_error("Fork failed.");
  }
  if (pid == 0) {
    // Child process
    close(pipefd[0]); // Close the read end of the pipe
    childProcess = true;
    signal(SIGALRM, handleAlarm);
    // Set an alarm
    alarm(static_cast<std::uint32_t>(childTimeout.count()));
    NASolver::Result result{};
    try {
      result = objective(arg);
    } catch (const std::exception& e) {
      alarm(0);
      close(pipefd[1]); // Close the write end of the pipe
      std::cerr << "[Optimizer]: Exception in objective function: " << e.what()
                << "\n";
      exit(1);
    }
    alarm(0); // Cancel the alarm
    // Write the result to the pipe
    const auto msg = result.yaml();
    write(pipefd[1], msg.c_str(), msg.size() + 1); // Write to the pipe
    close(pipefd[1]); // Close the write end of the pipe
    exit(0);          // End of child process
  }
  // Parent process
  close(pipefd[1]); // Close the write end of the pipe
  processData[pid] = {arg, pipefd[0]};
}

auto Optimizer::waitForChildProcess() -> void {
  int status = 0;
  // NOLINTNEXTLINE
  const pid_t pid = wait(&status);
  if (pid == -1) {
    killAllChildProcesses();
    throw std::runtime_error("Wait failed.");
  }
  if (WIFSIGNALED(status)) {
    if (!quiet) {
      std::cout << "k = " << processData[pid].arg << ": timeout\n";
    }
  } else if (WIFEXITED(status)) {
    // Read the result from the pipe
    std::stringstream msg;
    char              buffer[4096];
    ssize_t           bytesRead = 0;
    while ((bytesRead = read(processData[pid].readPipeFd, buffer,
                             sizeof(buffer) - 1)) > 0) {
      buffer[bytesRead] = '\0'; // Null-terminate the buffer
      msg << buffer;
    }
    if (bytesRead == -1) {
      killAllChildProcesses();
      throw std::runtime_error("Read failed.");
    }
    const auto& result = NASolver::Result::fromYAML(YAML::Load(msg.str()));
    if (result.isSat()) {
      if (!quiet) {
        std::cout << "k = " << processData[pid].arg << ": sat\n";
      }
      if (!minSat.has_value() || processData[pid].arg < minSat) {
        minSat   = processData[pid].arg;
        extremum = result;
      }
      for (const auto& [pPid, data] : processData) {
        if (data.arg > minSat) {
          if (!quiet) {
            std::cout << "k = " << data.arg << ": killing\n";
          }
          kill(pPid, SIGKILL);
        }
      }
    } else {
      if (!quiet) {
        std::cout << "k = " << processData[pid].arg << ": unsat\n";
      }
      if (!maxUnsat.has_value() || processData[pid].arg > maxUnsat) {
        maxUnsat = processData[pid].arg;
      }
      for (const auto& [pPid, data] : processData) {
        if (data.arg < maxUnsat) {
          if (!quiet) {
            std::cout << "k = " << data.arg << ": killing\n";
          }
          kill(pPid, SIGKILL);
        }
      }
    }
  } else {
    killAllChildProcesses();
    throw std::invalid_argument(
        "Exception occurred in child process, see message above.");
  }
  close(processData[pid].readPipeFd); // Close the read end of the pipe
  processData.erase(pid);
}

auto Optimizer::killAllChildProcesses() -> void {
  for (const auto& [pid, data] : processData) {
    kill(pid, SIGKILL);
    close(data.readPipeFd);
    waitpid(pid, nullptr, 0);
  }
  processData.clear();
}

auto Optimizer::minimize() -> void {
  const auto start = std::chrono::high_resolution_clock::now();
  // Initialize auxiliary variables
  extremum = std::nullopt;
  maxUnsat = std::nullopt;
  minSat   = std::nullopt;
  processData.clear();
  bool timeoutTriggered = false;
  for (std::size_t t = 10; !minSat.has_value(); t *= 10) {
    for (auto arg = maxUnsat.value_or(initialValue);
         !minSat.has_value() && arg <= maxValue; ++arg) {
      const auto delta = std::chrono::duration_cast<std::chrono::seconds>(
          timeout - (std::chrono::high_resolution_clock::now() - start));
      if (delta.count() <= 0) {
        timeoutTriggered = true;
        break; // only wait for running processes to terminate
      }
      forkChildProcess(arg, std::min(std::chrono::seconds(t), delta));
      if (getNSubProcsRunning() == maxNSubProcs) {
        waitForChildProcess();
      }
    }
  }
  while (isSubProcRunning()) {
    waitForChildProcess();
  }
  // no more child process is running
  assert(waitpid(-1, nullptr, WNOHANG) == -1 && errno == ECHILD);
  if (timeoutTriggered) {
    if (!quiet) {
      std::cerr
          << "[Optimizer]: Timeout happened while scanning for SAT instance.\n";
    }
    return;
  }
  if (!minSat.has_value()) {
    if (!quiet) {
      std::cerr << "[Optimizer]: No SAT instance found between initial value "
                   "and max value.\n";
    }
    return;
  }
  for (auto arg = minSat.value() - 1; maxUnsat.value_or(initialValue) < arg;
       --arg) {
    const auto delta = std::chrono::duration_cast<std::chrono::seconds>(
        timeout - (std::chrono::high_resolution_clock::now() - start));
    if (delta.count() <= 0) {
      timeoutTriggered = true;
      break; // only wait for running processes to terminate
    }
    assert(arg >= 0);
    forkChildProcess(static_cast<std::uint16_t>(arg), delta);
    if (getNSubProcsRunning() == maxNSubProcs) {
      waitForChildProcess();
    }
  }
  while (isSubProcRunning()) {
    waitForChildProcess();
  }
  // no more child process is running
  assert(waitpid(-1, nullptr, WNOHANG) == -1 && errno == ECHILD);
  if (timeoutTriggered && !quiet) {
    std::cerr
        << "[Optimizer]: Timeout happened while finding minimum SAT instance. "
           "Result may not be optimal.\n";
  }
}
} // namespace na
