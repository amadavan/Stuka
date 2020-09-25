//
// Created by Avinash Madavan on 2019-05-13.
//

#include <stuka/util/callback/progress.h>

stuka::util::callback::Progress::Progress(size_t n_max, const size_t precision, const bool carriage_return)
    : progress_bar_(n_max, precision, carriage_return) {}

void stuka::util::callback::Progress::callback(const stuka::OptimizeState state) {
  progress_bar_.update(state.nit);
}

void stuka::util::callback::Progress::finish(const stuka::OptimizeState state) {
  progress_bar_.close();
}
