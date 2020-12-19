//
// Created by Avinash Madavan on 1/15/19.
//

#include <stuka/dLP/benders_subproblem.h>

stuka::dLP::BendersSubproblem::BendersSubproblem(stuka::dLP::Subproblem &&sub, const Options &opts) : sub_(std::move(sub)), opts_(opts) {
  current_value_ = INF;

  LP::LinearProgram lp;
  lp.c = util::DenseOps::unique_copy(sub_.c);
  lp.A_ub = util::SparseOps::unique_copy(sub_.A_ub);
  lp.b_ub = util::DenseOps::unique_copy(sub_.b_ub);
  lp.A_eq = util::SparseOps::unique_copy(sub_.A_eq);
  lp.b_eq = util::DenseOps::unique_copy(sub_.b_eq);
  lp.lb = util::DenseOps::unique_copy(sub_.lb);
  lp.ub = util::DenseOps::unique_copy(sub_.ub);
  solver_ = util::createSolver(lp, opts);
}

stuka::dLP::BendersSubproblem::BendersSubproblem(BendersSubproblem &&sub) : opts_(sub.opts_),
                                                                            sub_{
                                                                                util::DenseOps::unique_copy(sub.sub_.c),
                                                                                util::SparseOps::unique_copy(sub.sub_.A_ub),
                                                                                util::DenseOps::unique_copy(sub.sub_.b_ub),
                                                                                util::SparseOps::unique_copy(sub.sub_.C_ub),
                                                                                util::SparseOps::unique_copy(sub.sub_.A_eq),
                                                                                util::DenseOps::unique_copy(sub.sub_.b_eq),
                                                                                util::SparseOps::unique_copy(sub.sub_.C_eq),
                                                                                util::DenseOps::unique_copy(sub.sub_.lb),
                                                                                util::DenseOps::unique_copy(sub.sub_.ub)} {
  solver_ = std::move(sub.solver_);
  current_value_ = sub.current_value_;
  
}

stuka::dLP::BendersCut stuka::dLP::BendersSubproblem::getBendersCut(const Eigen::VectorXd &x) {

  solver_->getLP().setRHS((sub_.b_ub) ? std::make_unique<Eigen::VectorXd>(*sub_.b_ub - *sub_.C_ub * x) : nullptr,
                          (sub_.b_eq) ? std::make_unique<Eigen::VectorXd>(*sub_.b_eq - *sub_.C_eq * x) : nullptr);

  // Solve the subproblem at the given point
  OptimizeState res;
  try {
    res = solver_->solve();
  } catch (std::exception &e) {
    throw std::runtime_error("getBendersCut: unable to solve subproblem");
  }
  if (res.status != 2) {
    std::cout << "Failed to solve subproblem: " << res.status << std::endl;
    throw std::runtime_error("getBendersCut: unable to solve subproblem");
  }

  Eigen::VectorXd aux = Eigen::VectorXd::Zero(x.size());
  if (sub_.b_ub)
    aux += sub_.C_ub->transpose() * res.dual_ub;
  if (sub_.b_eq)
    aux += sub_.C_eq->transpose() * res.dual_eq;

  current_value_ = res.fun;

  return {aux, -res.fun + aux.dot(x)};
}

const double &stuka::dLP::BendersSubproblem::getValue() {
  return current_value_;
}
