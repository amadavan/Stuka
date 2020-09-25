//
// Created by Avinash on 8/26/2020.
//

#ifndef STUKA_INCLUDE_STUKA_UTIL_PROGRESS_BAR_H_
#define STUKA_INCLUDE_STUKA_UTIL_PROGRESS_BAR_H_

#include <cstddef>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "timer.h"

namespace stuka {
namespace util {

template<typename T = size_t, typename = typename std::enable_if_t<std::is_arithmetic<T>::value, T>>
class ProgressBar {
 public:
  explicit ProgressBar(T n_max, const size_t precision = 3, const bool carriage_return = false)
      : n_max_(n_max), precision_(precision), carriage_return_(carriage_return) {
    size_t n = std::pow(10., precision_);
    n_p_ = ((size_t) (n_max_ / n) > 0) ? (n_max_ / n) : 1;

    print(0);
    timer_.start();
  }

  void update(T index) {
    if (index % n_p_ == 0)
      print(100 * (double) index / n_max_);
  }

  void close() {
    print(100.);
    if (carriage_return_) std::cout << std::endl;
  }

 private:
  const T n_max_;
  T n_p_;
  const size_t precision_;
  const bool carriage_return_;

  Timer timer_;

  void print(double percentage) {
    double elapsed = timer_.elapsed().count();
    double estimated = (100. - percentage) * (elapsed / percentage);

    std::cout << std::right << std::setw(precision_ + 2) << std::setprecision(precision_) << percentage << "%"
              << "          "
              << "Elapsed: " << std::left << std::setw(10) << std::setprecision(10) << elapsed << "s"
              << "          "
              << "Estimated: " << std::left << std::setw(10) << std::setprecision(10) << estimated << "s"
        ;
    if (carriage_return_) std::cout << '\r' << std::flush;
    else std::cout << std::endl;
  }
};

}
}

#endif //STUKA_INCLUDE_STUKA_UTIL_PROGRESS_BAR_H_
