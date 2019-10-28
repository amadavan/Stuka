//
// Created by Avinash Madavan on 12/4/18.
//

#ifndef STUKA_LP_GUROBI_H
#define STUKA_LP_GUROBI_H

#include <gurobi_c++.h>

#include "base_solver.h"
#include "lp.h"
#include "gurobi_lp.h"

namespace stuka { namespace LP {

class GurobiSolver : public BaseLPSolver {
 public:
  explicit GurobiSolver(const LinearProgram &lp, const Options &opts = Options());

  ~GurobiSolver() override = default;

  BaseLinearProgram &getLP() override;

  void iterate() override;

  bool terminate() override;

  const stuka::OptimizeState getState() override;

 private:
  GurobiLinearProgram prog_;
};

}}

#endif //STUKA_LP_GUROBI_H
