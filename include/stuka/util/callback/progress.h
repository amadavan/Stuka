//
// Created by Avinash Madavan on 2019-05-13.
//

#ifndef STUKA_CALLBACK_PROGRESS_H
#define STUKA_CALLBACK_PROGRESS_H

#include <iostream>
#include <iomanip>

#include "base_callback.h"
#include "../timer.h"

namespace stuka { namespace util { namespace callback {
class Progress : public BaseCallback {
 public:
  Progress(size_t n_max_);

  void callback(const OptimizeState state) override;

  void finish(const OptimizeState state) override;

 private:
  size_t n_max_;
  size_t n_p_;

  Timer timer_;
};
}}}

#endif //STUKA_CALLBACK_PROGRESS_H
