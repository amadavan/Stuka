//
// Created by Avinash Madavan on 12/1/18.
//

#ifndef STUKA_LP_BASE_SOLVER_H
#define STUKA_LP_BASE_SOLVER_H

#ifndef DEFAULT_MAX_ITER
#define DEFAULT_MAX_ITER 100
#endif

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "../base_solver.h"
#include "base_lp.h"

namespace stuka { namespace LP {

  class BaseLPSolver : public BaseSolver {
  public:
    BaseLPSolver(const LinearProgram &lp, const Options &opts) : BaseSolver(opts) {}

    virtual ~BaseLPSolver() override = default;

    /* Get the BaseLP of the solver.
     *
     * Each solver may use a reformulated version of the original LP that is given in the constructor. This function
     * returns the reformulated LP as a BaseLP to allow for updates to the solver. This step also provides a method to
     * trigger any updates that may be necessary due to changes in the LP.
     */
    virtual BaseLinearProgram &getLP() = 0;
  };

}}

#endif //STUKA_LP_BASE_SOLVER_H
