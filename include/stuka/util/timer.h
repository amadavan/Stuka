//
// Created by Avinash Madavan on 1/21/19.
//

#ifndef STUKA_UTIL_TIMER_H
#define STUKA_UTIL_TIMER_H

#include <chrono>

namespace stuka { namespace util {
class Timer {
 public:
  Timer();

  ~Timer() = default;

  void start();

  std::chrono::duration<double> elapsed();

 protected:
  std::chrono::time_point<std::chrono::high_resolution_clock> time_;
};
}}

#endif //STUKA_UTIL_TIMER_H
