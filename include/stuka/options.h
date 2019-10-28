//
// Created by Avinash Madavan on 1/15/19.
//

#ifndef STUKA_OPTIONS_H
#define STUKA_OPTIONS_H

#include <stddef.h>

#include <functional>
#include <memory>

#include <Eigen/Core>

#include "constants.h"
#include "util/callback/base_callback.h"

namespace stuka {
struct Options {
  // General options
  Eigen::VectorXd x0;                                           // Initial state value
  double tol = 1e-8;                                            // Solver tolerance
  size_t max_iter = 0;                                          // Maximum number of iterations
  double step = 1.;                                             // Step size for gradient descent type algorithms

  // For linear and quadratic programs
  Eigen::VectorXd dual_ub0;                                     // Initial upper-bound dual variable value
  Eigen::VectorXd dual_eq0;                                     // Initial equality-constraint dual variable value

  // For decomposed linear programs
  stuka::Solver lp_solver = DEFAULT_LP_SOLVER;                  // Solver to use for linear programs
  stuka::Solver qp_solver = DEFAULT_QP_SOLVER;                  // Solver to use for quadratic programs
  stuka::Solver dlp_solver = DEFAULT_DLP_SOLVER;                // Solver to use for decomposed linear programs
  stuka::Solver stochastic_solver = DEFAULT_STOCHASTIC_SOLVER;  // Solver to use for stochastic programs

  // Callback
  std::shared_ptr<util::callback::BaseCallback> callback = nullptr;

  // Hack for cleaner python implementation
  virtual std::string name() { return "Options"; }
};
}

#endif //STUKA_OPTIONS_H
