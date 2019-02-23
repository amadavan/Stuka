//
// Created by Avinash Madavan on 1/7/19.
//

#ifndef STUKA_QP_BASE_SOLVER_H
#define STUKA_QP_BASE_SOLVER_H

#ifndef DEFAULT_MAX_ITER
#define DEFAULT_MAX_ITER 100
#endif

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "../base_solver.h"
#include "base_qp.h"

namespace stuka { namespace QP {

  class BaseQPSolver : public BaseSolver {
  public:
    BaseQPSolver(const QuadraticProgram &prog, const Options &opts) : BaseSolver(opts) {}

    virtual ~BaseQPSolver() override = default;

    /* Get the BaseLP of the solver.
     *
     * Each solver may use a reformulated version of the original LP that is given in the constructor. This function
     * returns the reformulated LP as a BaseLP to allow for updates to the solver. This step also provides a method to
     * trigger any updates that may be necessary due to changes in the LP.
     */
    virtual BaseQuadraticProgram &getQP() = 0;
  };

}}

#endif //STUKA_QP_BASE_SOLVER_H
