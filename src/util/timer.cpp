//
// Created by Avinash Madavan on 1/21/19.
//

#include <stuka/util/timer.h>

stuka::util::Timer::Timer() {
  time_ = std::chrono::high_resolution_clock::now();
}

void stuka::util::Timer::start() {
  time_ = std::chrono::high_resolution_clock::now();
}

std::chrono::duration<double> stuka::util::Timer::elapsed() {
  std::chrono::time_point<std::chrono::high_resolution_clock> ctime = std::chrono::high_resolution_clock::now();
  return {ctime - time_};
}
