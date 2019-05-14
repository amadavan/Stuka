//
// Created by Avinash Madavan on 2019-05-13.
//

#ifndef STUKA_CALLBACK_COMPOSITE_H
#define STUKA_CALLBACK_COMPOSITE_H

#include <memory>
#include <vector>

#include "base_callback.h"

namespace stuka { namespace util { namespace callback {
  class Composite : public BaseCallback {
  public:
    Composite(const std::vector<std::shared_ptr<BaseCallback>> &cbs_);

    void initialize(const OptimizeState state) override;

    void callback(const OptimizeState state) override;

  public:

  private:
    const std::vector<std::shared_ptr<BaseCallback>> cbs_;
  };
}}}

#endif //STUKA_CALLBACK_COMPOSITE_H
