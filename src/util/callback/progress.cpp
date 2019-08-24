//
// Created by Avinash Madavan on 2019-05-13.
//

#include <stuka/util/callback/progress.h>

stuka::util::callback::Progress::Progress(size_t n_max) : n_max_(n_max) {
  n_p_ = ((size_t) (n_max_ / 10000) > 0) ? (n_max_ / 10000) : 1;
  timer_.start();
}

void stuka::util::callback::Progress::callback(const stuka::OptimizeState state) {
  if (state.nit % n_p_ == 0) {
    double cp = state.nit * (100. / n_max_);
    double elapsed = timer_.elapsed().count();
    std::cout << std::right << std::setw(5) << std::setprecision(4) << cp << "%"
              << "          "
              << "Elapsed: " << std::left << std::setw(10) << std::setprecision(10) << elapsed << "s"
              << "          "
              << "Estimated: " << std::left << std::setw(10) << std::setprecision(10) << (100. - cp) * (elapsed / cp)
              << "s"
              << std::endl;
  }
}

void stuka::util::callback::Progress::finish(const stuka::OptimizeState state) {
  double elapsed = timer_.elapsed().count();
  std::cout << std::right << std::setw(5) << std::setprecision(4) << 100 << "%"
            << "          "
            << "Total Time: " << std::left << std::setw(10) << std::setprecision(10) << elapsed << "s"
            << std::endl;
}
