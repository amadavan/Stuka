//
// Created by Avinash Madavan on 2019-05-13.
//

#ifndef STUKA_CALLBACK_PROGRESS_H
#define STUKA_CALLBACK_PROGRESS_H

#include <iostream>

#include "base_callback.h"

namespace stuka { namespace util { namespace callback {
  class Progress : public BaseCallback {
  public:
    Progress();

    Progress(unsigned int n_max_);

    void callback(const OptimizeState state) override;

  private:
    unsigned int n_max_;
    unsigned int n_p_;
    double p_;
  };
}}}

#endif //STUKA_CALLBACK_PROGRESS_H
