//
// Created by Avinash on 8/31/2020.
//

#ifndef STUKA_INCLUDE_STUKA_DLP_NAIVE_BENDERS_H_
#define STUKA_INCLUDE_STUKA_DLP_NAIVE_BENDERS_H_

#include "../constants.h"
#include "../util/solver_factory.h"
#include "../LP/base_solver.h"
#include "../util/dense_ops.h"
#include "../util/sparse_ops.h"

#include "base_solver.h"
#include "benders_cut.h"
#include "naive_benders_subproblem.h"

namespace stuka { namespace dLP {
class NaiveBendersDecomposition : public BaseDLPSolver {
 public:
  explicit NaiveBendersDecomposition(const DecomposedLinearProgram &dlp, const Options &opts = Options());

  ~NaiveBendersDecomposition() override = default;

  void iterate() override;

  bool terminate() override;

  const OptimizeState getState() override;

 private:
  OptimizeState solveMasterProblem(const std::vector<BendersCut> &cuts);

  Eigen::VectorXd x_;                                     ///< Current iterate
  size_t n_sub_;
  size_t n_dim_master_;
  size_t n_sub_calls_;

  const Options opts_;                                    ///< Options

  std::vector<NaiveBendersSubproblem> subproblems_;       ///< Subproblems
  Eigen::VectorXd subproblem_values_;                     ///< Function value of each subproblem

  LP::LinearProgram master_lp_;
};

}}
#endif //STUKA_INCLUDE_STUKA_DLP_NAIVE_BENDERS_H_
