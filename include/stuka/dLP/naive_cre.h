//
// Created by Avinash on 8/31/2020.
//

#ifndef STUKA_INCLUDE_STUKA_DLP_NAIVE_CRE_H_
#define STUKA_INCLUDE_STUKA_DLP_NAIVE_CRE_H_

#include <mutex>          // std::mutex

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>

#include "../options.h"
#include "../constants.h"
#include "../util/solver_factory.h"
#include "../util/dense_ops.h"
#include "../util/sparse_ops.h"
#include "../LP/base_solver.h"
#include "../QP/base_solver.h"

#include "base_solver.h"
#include "naive_cre_subproblem.h"
#include "critical_region.h"

namespace stuka { namespace dLP {
class NaiveCRE : public BaseDLPSolver {
 public:
  explicit NaiveCRE(const DecomposedLinearProgram &dlp, const Options &opts = Options());

  ~NaiveCRE() override = default;

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

  Eigen::SparseMatrix<double> A_activeT_;
  Eigen::MatrixXd previous_directions_;

  LP::LinearProgram master_lp_;                         ///< Master problem parameters
  QP::QuadraticProgram residual_qp_;                    ///< Residual problem
//  std::unique_ptr<LP::BaseLPSolver> master_solver_;             // Solver for master problem
//  std::unique_ptr<QP::BaseQPSolver> residual_solver_;           // Solver for determining residual
  // TODO: remove residual_solver

  std::vector<NaiveCRESubproblem> subproblems_;

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
#endif //STUKA_INCLUDE_STUKA_DLP_NAIVE_CRE_H_
