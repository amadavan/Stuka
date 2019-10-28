//
// Created by Avinash Madavan on 2019-06-26.
//

#ifndef STUKA_LP_BASE_PRECONDITIONER_H
#define STUKA_LP_BASE_PRECONDITIONER_H

#include "../base_lp.h"

namespace stuka { namespace LP {
class BasePreconditioner : public BaseLinearProgram {
 public:
  BasePreconditioner(const std::shared_ptr<BaseLinearProgram> &next) : next_(next) {
    if (!next_) throw std::runtime_error("Preconditioner requires program.");
  }

  const std::shared_ptr<BaseLinearProgram> &next() { return next_; }

 private:
  const std::shared_ptr<BaseLinearProgram> next_;
};
}}

#endif //STUKA_LP_BASE_PRECONDITIONER_H
