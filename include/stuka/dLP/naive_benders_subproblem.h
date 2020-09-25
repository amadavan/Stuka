//
// Created by Avinash on 8/31/2020.
//

#ifndef STUKA_INCLUDE_STUKA_DLP_NAIVE_BENDERS_SUBPROBLEM_H_
#define STUKA_INCLUDE_STUKA_DLP_NAIVE_BENDERS_SUBPROBLEM_H_

#include <memory>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "../constants.h"
#include "../util/solver_factory.h"
#include "../LP/base_solver.h"
#include "subproblem.h"
#include "benders_cut.h"

namespace stuka { namespace dLP {
class NaiveBendersSubproblem {
 public:
  explicit NaiveBendersSubproblem(Subproblem sub, const Options = Options());
  NaiveBendersSubproblem(NaiveBendersSubproblem &&sub) noexcept;

  BendersCut getBendersCut(const Eigen::VectorXd &x);

  const double &getValue();
 private:
  std::unique_ptr<LP::BaseLPSolver> solver_;
  LP::LinearProgram prog_;
  const Options opts_;
  const Subproblem sub_;
  double current_value_;
};
}}

#endif //STUKA_INCLUDE_STUKA_DLP_NAIVE_BENDERS_SUBPROBLEM_H_
