//
// Created by Avinash Madavan on 1/10/19.
//

#ifndef STUKA_OPTIMIZE_STATE_H
#define STUKA_OPTIMIZE_STATE_H

#include <Eigen/Core>

namespace stuka {

struct OptimizeState {
  Eigen::VectorXd x;            ///< Primal variable
  Eigen::VectorXd dual_ub;      ///< Dual variable associated with inequality constraints
  Eigen::VectorXd dual_eq;      ///< Dual variable associated with equality constraints
  Eigen::VectorXd dual_x_lb;    ///< Dual variable associated with primal variable lower bounds
  Eigen::VectorXd dual_x_ub;    ///< Dual variable associated with primal variable upper bounds

  double fun = 0;         ///< Function value
  double error = 0;       ///< Solver error
  int status = 0;         ///< Status code
  size_t nit = 0;         ///< Iteration count
  size_t nit_sub = 0;     ///< Number of computations of subproblems for decomposed linear programs
  size_t nit_con_gen = 0; ///< Number of constraint generation iterations

  double runtime = 0;     ///< Runtime (in s)
};

}

#endif //STUKA_OPTIMIZE_STATE_H
