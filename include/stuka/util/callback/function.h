//
// Created by Avinash Madavan on 2019-05-08.
//

#ifndef STUKA_CALLBACK_FUNCTION_H
#define STUKA_CALLBACK_FUNCTION_H

#include "base_callback.h"

namespace stuka { namespace util { namespace callback {
  class Function : public BaseCallback {
  public:
    Function(const std::function<void(const OptimizeState)> func);

    void callback(const OptimizeState state) override;

  private:
    const std::function<void(const OptimizeState state)> func_;
  };
}}}

#endif //STUKA_CALLBACK_FUNCTION_H
