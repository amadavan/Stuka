//
// Created by Avinash Madavan on 1/10/19.
//

#ifndef STUKA_dLP_BASE_SOLVER_H
#define STUKA_dLP_BASE_SOLVER_H

#ifndef DEFAULT_MAX_ITER
#define DEFAULT_MAX_ITER 100
#endif

#include "../base_solver.h"
#include "decomposed_lp.h"

namespace stuka { namespace dLP {
class BaseDLPSolver : public BaseSolver {
 public:
  BaseDLPSolver(const DecomposedLinearProgram &dlp, const Options &opts) : BaseSolver(opts) {}

  virtual ~BaseDLPSolver() override = default;
};
}}

#endif //STUKA_dLP_BASE_SOLVER_H
