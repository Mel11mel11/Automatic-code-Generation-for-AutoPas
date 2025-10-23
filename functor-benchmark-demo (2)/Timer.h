#pragma once

#include <chrono>
#include <iomanip>
#include <sstream>
#include <string>

using namespace std::chrono;

/**
 * Timer class to stop times.
 */
class Timer {
 public:
  Timer() = default;

  virtual ~Timer() = default;

  /**
   * start the timer.
   */
  void start() {
    if (_currentlyRunning) {
        throw("Trying to start a timer that is already started!");
    }
    _currentlyRunning = true;
    _startTime = high_resolution_clock::now();
  }

  /**
   * Stops the timer and returns the time elapsed in nanoseconds since the last call to start.
   * It also adds the duration to the total time.
   * @return elapsed time in nanoseconds
   */
  long stop() {
    const auto time(high_resolution_clock::now());

    if (not _currentlyRunning) {
        throw("Trying to stop a timer that was not started!");
    }
    _currentlyRunning = false;

    const auto diff = duration_cast<nanoseconds>(time - _startTime).count();

    _totalTime += diff;
    return diff;
  }

  /**
   * Resets the timer to 0.
   */
  void reset() {_totalTime = 0;}

  /**
   * Adds the given amount of nanoseconds to the total time.
   * @param nanoseconds
   */
  void addTime(long nanoseconds) {_totalTime += nanoseconds;}

  /**
   * Get total accumulated time.
   * @return Total time in nano seconds.
   */
  [[nodiscard]] long getTotalTime() const { return _totalTime; }

  /**
   * Create a date stamp for the current moment with the given format.
   * @param format Date stamp format.
   * @return String representation of the current date in the given format.
   */
  static std::string getDateStamp(const std::string &format = "%Y-%m-%d_%H-%M-%S") {
    auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::ostringstream nowStrStr;
    tm unused;
    std::stringstream ss;
    ss << std::put_time(localtime_r(&now, &unused), format.c_str());
    return ss.str();
  }

 private:
  /**
   * Time point of last call of start().
   */
  std::chrono::high_resolution_clock::time_point _startTime{};

  /**
   * Accumulated total time.
   */
  long _totalTime = 0;

  /**
   * Indicator if this timer currently is measuring.
   */
  bool _currentlyRunning = false;
};
