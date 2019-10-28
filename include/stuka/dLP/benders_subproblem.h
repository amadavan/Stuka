//
// Created by Avinash Madavan on 1/15/19.
//

#ifndef STUKA_dLP_BENDERS_SUBPROBLEM_H
#define STUKA_dLP_BENDERS_SUBPROBLEM_H

#include <memory>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "../constants.h"
#include "stuka/util/solver_factory.h"
#include "../LP/base_solver.h"
#include "subproblem.h"
#include "benders_cut.h"

namespace stuka { namespace dLP {
class BendersSubproblem {
 public:
  explicit BendersSubproblem(Subproblem sub, const Options = Options());
  BendersSubproblem(BendersSubproblem &&sub) noexcept;

  BendersCut getBendersCut(const Eigen::VectorXd &x);

  const double &getValue();
 private:
  std::unique_ptr<LP::BaseLPSolver> solver_;
  const Subproblem sub_;
  double current_value_;
};
}}

#endif //STUKA_dLP_BENDERS_SUBPROBLEM_H
