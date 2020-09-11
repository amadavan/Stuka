//
// Created by Avinash Madavan on 1/15/19.
//

#ifndef STUKA_dLP_BENDERS_H
#define STUKA_dLP_BENDERS_H

#include "../constants.h"
#include "../util/solver_factory.h"
#include "../LP/base_solver.h"

#include "base_solver.h"
#include "benders_cut.h"
#include "benders_subproblem.h"

namespace stuka { namespace dLP {
class BendersDecomposition : public BaseDLPSolver {
 public:
  explicit BendersDecomposition(const DecomposedLinearProgram &dlp, const Options &opts = Options());

  ~BendersDecomposition() override = default;

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

  std::unique_ptr<LP::BaseLPSolver> master_solver_;       ///< Solver for master problem
  std::vector<BendersSubproblem> subproblems_;            ///< Subproblems
  Eigen::VectorXd subproblem_values_;                     ///< Function value of each subproblem

  LP::LinearProgram master_lp_;
};

}}

#endif //STUKA_dLP_BENDERS_H
