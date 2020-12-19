//
// Created by Avinash on 8/19/2020.
//

#ifndef STUKA_INCLUDE_STUKA_LP_LAZY_CONSTRAINTS_H_
#define STUKA_INCLUDE_STUKA_LP_LAZY_CONSTRAINTS_H_

#include "base_solver.h"
#include "lazy_lp.h"
#include "../util/solver_factory.h"

namespace stuka { namespace LP {

class LazyConstraints : public BaseLPSolver {
 public:
  LazyConstraints(const LinearProgram &lp, const Options &opts);
  ~LazyConstraints() override = default;

  void iterate() override;

  bool terminate() override;

  const OptimizeState getState() override;

  BaseLinearProgram &getLP() override;

 private:
  size_t nit_;
  size_t constraints_added_;

  std::unique_ptr<BaseLPSolver> solver_;
  LazyLinearProgram lp_;

  OptimizeState state_;
};

}}

#endif //STUKA_INCLUDE_STUKA_LP_LAZY_CONSTRAINTS_H_
