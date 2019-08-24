//
// Created by Avinash Madavan on 2019-05-07.
//

#ifndef STUKA_STOCHASTIC_BASE_SOLVER_H
#define STUKA_STOCHASTIC_BASE_SOLVER_H

#ifndef DEFAULT_MAX_ITER
#define DEFAULT_MAX_ITER 100
#endif

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "../base_solver.h"
#include "program.h"

namespace stuka { namespace stochastic {

  class BaseStochasticSolver : public BaseSolver {
  public:
    BaseStochasticSolver(const Program &prog, const Options &opts) : BaseSolver(opts) {}

    virtual ~BaseStochasticSolver() override = default;
  };

}}

#endif //STUKA_STOCHASTIC_BASE_SOLVER_H
