//
// Created by Avinash Madavan on 1/7/19.
//

#ifndef STUKA_QP_GUROBI_H
#define STUKA_QP_GUROBI_H

#include <gurobi_c++.h>

#include "../optimize_state.h"

#include "base_solver.h"
#include "qp.h"
#include "gurobi_qp.h"

namespace stuka { namespace QP {

  class GurobiSolver : public BaseQPSolver {
  public:
    explicit GurobiSolver(const QuadraticProgram &qp, const Options &opts = Options());

    ~GurobiSolver() override = default;

    BaseQuadraticProgram &getQP() override;

    void iterate() override;

    bool terminate() override;

    const stuka::OptimizeState getState() override;

  private:
    GurobiQuadraticProgram prog_;
  };

}}

#endif //STUKA_QP_GUROBI_H
