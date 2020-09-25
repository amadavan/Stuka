//
// Created by Avinash Madavan on 2019-05-13.
//

#ifndef STUKA_CALLBACK_PROGRESS_H
#define STUKA_CALLBACK_PROGRESS_H

#include <iostream>
#include <iomanip>

#include "base_callback.h"
#include "../timer.h"
#include "../progress_bar.h"

namespace stuka { namespace util { namespace callback {
class Progress : public BaseCallback {
 public:
  Progress(size_t n_max_, const size_t precision = 3, const bool carriage_return = false);

  void callback(const OptimizeState state) override;

  void finish(const OptimizeState state) override;

 private:
  ProgressBar<size_t> progress_bar_;
};
}}}

#endif //STUKA_CALLBACK_PROGRESS_H
