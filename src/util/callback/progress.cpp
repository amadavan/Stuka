//
// Created by Avinash Madavan on 2019-05-13.
//

#include <stuka/util/callback/progress.h>

stuka::util::callback::Progress::Progress(unsigned int n_max) : n_max_(n_max) {
  n_p_ = ((int) (n_max_ / 100) > 0) ? (n_max_ / 100) : 1;
}

void stuka::util::callback::Progress::callback(const stuka::OptimizeState state) {
  std::cout << state.nit << std::endl;
  std::cout << n_p_ << std::endl;
  if (state.nit % n_p_ == 0) {
    std::cout << (int) (state.nit * 100 / n_max_) << "%" << std::endl;
  }
}
