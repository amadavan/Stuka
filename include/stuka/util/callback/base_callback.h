//
// Created by Avinash Madavan on 2019-05-08.
//

#ifndef STUKA_CALLBACK_BASE_CALLBACK_H
#define STUKA_CALLBACK_BASE_CALLBACK_H

#include "stuka/optimize_state.h"

namespace stuka { namespace util { namespace callback {
  class BaseCallback {
  public:
    virtual ~BaseCallback() = default;

    virtual void initialize(const OptimizeState state) {};

    virtual void callback(const OptimizeState state) = 0;

    virtual void finish(const OptimizeState state) {};
  };
}}}

#endif //STUKA_CALLBACK_BASE_CALLBACK_H
