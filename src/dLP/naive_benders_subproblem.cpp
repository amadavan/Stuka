//
// Created by Avinash on 8/31/2020.
//

//
// Created by Avinash Madavan on 1/15/19.
//

#include <stuka/dLP/naive_benders_subproblem.h>

stuka::dLP::NaiveBendersSubproblem::NaiveBendersSubproblem(stuka::dLP::Subproblem &&sub, const Options &opts) : sub_(std::move(sub)), opts_(opts) {
  current_value_ = INF;

  prog_.c = util::DenseOps::unique_copy(sub_.c);
  prog_.A_ub = util::SparseOps::unique_copy(sub_.A_ub);
  prog_.b_ub = util::DenseOps::unique_copy(sub_.b_ub);
  prog_.A_eq = util::SparseOps::unique_copy(sub_.A_eq);
  prog_.b_eq = util::DenseOps::unique_copy(sub_.b_eq);
  prog_.lb = util::DenseOps::unique_copy(sub_.lb);
  prog_.ub = util::DenseOps::unique_copy(sub_.ub);
}

stuka::dLP::NaiveBendersSubproblem::NaiveBendersSubproblem(NaiveBendersSubproblem &&sub) : opts_(sub.opts_),
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
  current_value_ = sub.current_value_;

  prog_.c = util::DenseOps::unique_copy(sub.sub_.c);
  prog_.A_ub = util::SparseOps::unique_copy(sub.sub_.A_ub);
  prog_.b_ub = util::DenseOps::unique_copy(sub.sub_.b_ub);
  prog_.A_eq = util::SparseOps::unique_copy(sub.sub_.A_eq);
  prog_.b_eq = util::DenseOps::unique_copy(sub.sub_.b_eq);
  prog_.lb = util::DenseOps::unique_copy(sub.sub_.lb);
  prog_.ub = util::DenseOps::unique_copy(sub.sub_.ub);
  
}

stuka::dLP::BendersCut stuka::dLP::NaiveBendersSubproblem::getBendersCut(const Eigen::VectorXd &x) {

  if (sub_.b_ub) prog_.b_ub = std::make_unique<Eigen::VectorXd>(*sub_.b_ub - *sub_.C_ub * x);
  if (sub_.b_eq) prog_.b_eq = std::make_unique<Eigen::VectorXd>(*sub_.b_eq - *sub_.C_eq * x);

  std::unique_ptr<LP::BaseLPSolver> solver = stuka::util::createSolver(prog_, opts_);

  // Solve the subproblem at the given point
  OptimizeState res;
  try {
    res = solver->solve();
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

const double &stuka::dLP::NaiveBendersSubproblem::getValue() {
  return current_value_;
}
