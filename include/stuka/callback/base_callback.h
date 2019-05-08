//
// Created by Avinash Madavan on 2019-05-08.
//

#ifndef STUKA_CALLBACK_BASE_CALLBACK_H
#define STUKA_CALLBACK_BASE_CALLBACK_H

#include "../optimize_state.h"

namespace stuka { namespace callback {
  class BaseCallback {
  public:
    virtual ~BaseCallback() = default;

    virtual void callback(OptimizeState state) = 0;
  };
}}

#endif //STUKA_CALLBACK_BASE_CALLBACK_H
