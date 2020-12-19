//
// Created by Avinash Madavan on 1/15/19.
//

#ifndef STUKA_OPTIONS_H
#define STUKA_OPTIONS_H

#include <stddef.h>
#include <functional>
#include <memory>

#include <Eigen/Core>
#include <gurobi_c++.h>

#include "constants.h"
#include "defaults.h"
#include "util/callback/base_callback.h"
#include "stochastic/measure/measure.h"
#include "stochastic/measure/expected_value.h"

namespace stuka {
/** Structure for all user options.
 *
 * Contains all necessary information for solvers to initialize their routines and perform the optimization. This
 * includes methods to select the optimization routine, as well as
 */
struct Options {
  // General options
  Eigen::VectorXd x0;                                           ///< Initial state value
  double tol = 1e-8;                                            ///< Solver tolerance
  size_t max_iter = 0;                                          ///< Maximum number of iterations
  double step = 1.;                                             ///< Step size for gradient descent type algorithms

  // For linear and quadratic programs
  Eigen::VectorXd dual_ub0;                                     ///< Initial upper-bound dual variable value
  Eigen::VectorXd dual_eq0;                                     ///< Initial equality-constraint dual variable value

  // For decomposed linear programs

  // For stochastic programs
  /// Measure to apply for stochastic optimization
  std::shared_ptr<stuka::stochastic::ProgramMeasure> measure = std::make_shared<stuka::stochastic::ExpectedValue>();

  // Solver specification
  stuka::Solver lp_solver = DEFAULT_LP_SOLVER;                  ///< Solver to use for linear programs
  stuka::Solver qp_solver = DEFAULT_QP_SOLVER;                  ///< Solver to use for quadratic programs
  stuka::Solver dlp_solver = DEFAULT_DLP_SOLVER;                ///< Solver to use for decomposed linear programs
  stuka::Solver stochastic_solver = DEFAULT_STOCHASTIC_SOLVER;  ///< Solver to use for stochastic programs

  // Other solver specific options
  double cre_step = DEFAULT_CRE_STEP;                           ///< Step size for CRE algorithm

  // Lazy constraint generation (CURRENTLY ONLY FOR LINEAR PROGRAMS)
  bool lazy = false;
  std::vector<size_t> active_indices = std::vector<size_t>();   ///< Indices of UB constraints that begin active

  // Subproblem options
  bool lazy_sub = false;

  // Gurobi solver preferences
  GurobiMethod gurobi_method = GurobiMethod::CONCURRENT;
  // todo: implement for qp

  // Callback
  std::shared_ptr<util::callback::BaseCallback> callback = nullptr;

 private:
  // Hack for cleaner python implementation
  virtual std::string name() { return "Options"; }
};
}

#endif //STUKA_OPTIONS_H
