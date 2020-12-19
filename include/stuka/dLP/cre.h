//
// Created by Avinash Madavan on 1/21/19.
//

#ifndef STUKA_dLP_CRE_H
#define STUKA_dLP_CRE_H

#include <mutex>          // std::mutex

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>

#include "../options.h"
#include "../util/solver_factory.h"
#include "../util/dense_ops.h"
#include "../util/sparse_ops.h"
#include "../LP/base_solver.h"
#include "../QP/base_solver.h"

#include "base_solver.h"
#include "cre_subproblem.h"
#include "critical_region.h"

namespace stuka { namespace dLP {
class CRE : public BaseDLPSolver {
 public:
  explicit CRE(const DecomposedLinearProgram &dlp, const Options &opts = Options());

  ~CRE() override = default;

  void iterate() override;

  bool terminate() override;

  const OptimizeState getState() override;

 private:
  std::mutex mtx_cre_;
  Eigen::VectorXd x_;                 ///< Current iterate
  size_t n_sub_;
  long n_dim_master_;
  size_t n_active_;
  size_t n_visit_;
  size_t n_sub_calls_;
  bool first_run_;

  const Options opts_;                ///< Options

  Eigen::VectorXd a_opt_;             ///< Residual vector
  double opt_val_;                    ///< Current optimal value
  Eigen::VectorXd x_opt_;             ///< Current optimum

  CriticalRegion cr_;                 ///< Global critical region
  long n_ub_cr_;

  // Master problem properties
  size_t n_eq_;
  size_t n_ub_;
  Eigen::Array<bool, Eigen::Dynamic, 1> bone_eq_;
  Eigen::SparseMatrix<double> eye_;
  Eigen::SparseMatrix<double> negative_eye_;

  LP::LinearProgram master_lp_;                                 // Master problem parameters
  std::unique_ptr<LP::BaseLPSolver> master_solver_;             // Solver for master problem
  std::unique_ptr<QP::BaseQPSolver> residual_solver_;           // Solver for determining residual

  std::vector<CRESubproblem> subproblems_;

  struct CriticalRegionData {
    CriticalRegion cr;
    size_t ind_constr = 0;
    size_t ind_add = 0;
  };

  std::vector<CriticalRegionData> critical_regions_;

  OptimizeState solveMasterProblem();
  const Eigen::VectorXd computeResidualFirst(const OptimizeState &x);
  const Eigen::VectorXd computeResidualSubsequent();
  void removeCriticalRegion(const CriticalRegionData &crdata);
  void addCriticalRegion(CriticalRegionData &crdata);
};
}}

#endif //STUKA_dLP_CRE_H
