//
// Created by Avinash Madavan on 1/15/19.
//

#ifndef STUKA_dLP_CRE_SUBPROBLEM_H
#define STUKA_dLP_CRE_SUBPROBLEM_H

#include <memory>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SPQRSupport>

#include "../util/solver_factory.h"
#include "../LP/base_solver.h"
#include "../LP/lp.h"
#include "subproblem.h"
#include "critical_region.h"

namespace stuka { namespace dLP {
class CRESubproblem {
 public:
  explicit CRESubproblem(Subproblem sub, const Options = Options());

  CriticalRegion computeCriticalRegion(const Eigen::VectorXd &x) const;
 private:
  std::unique_ptr<LP::BaseLPSolver> solver_;

  long n_dim_;
  long n_dim_master_;
  long n_bounds_;

  const Subproblem sub_;
};
}}

#endif //STUKA_dLP_CRE_SUBPROBLEM_H
