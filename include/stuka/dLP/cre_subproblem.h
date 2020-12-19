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

#include "../options.h"
#include "../constants.h"
#include "../util/sparse_ops.h"
#include "../util/dense_ops.h"
#include "../util/solver_factory.h"
#include "../LP/base_solver.h"
#include "../LP/lp.h"
#include "subproblem.h"
#include "critical_region.h"

namespace stuka { namespace dLP {
class CRESubproblem {
 public:
  explicit CRESubproblem(Subproblem &&sub, const Options & = Options());
  CRESubproblem(CRESubproblem &&sub);

  [[nodiscard]] CriticalRegion computeCriticalRegion(const Eigen::VectorXd &x);
 private:
  std::unique_ptr<LP::BaseLPSolver> solver_;

  long n_dim_;
  long n_dim_master_;
  long n_bounds_;

  Eigen::VectorXd b_lb_;
  Eigen::VectorXd b_ub_;
  util::DenseOps::ActiveSet bone_eq_;           // Array of ones whose length is the number of equality constraints
  util::DenseOps::ActiveSet xlb_constraints_;   // Boolean array denoting whether lb constraint exists
  util::DenseOps::ActiveSet xub_constraints_;   // Boolean array denoting whether ub constraint exists
  Eigen::SparseMatrix<double> A_xlb_;           // Sparse identity matrix corresponding to lb constraints
  Eigen::SparseMatrix<double> A_xub_;           // Sparse identity matrix corresponding to ub constraints

  const Subproblem sub_;
  const Options opts_;
};
}}

#endif //STUKA_dLP_CRE_SUBPROBLEM_H
