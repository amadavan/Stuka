//
// Created by Avinash on 8/31/2020.
//

//
// Created by Avinash Madavan on 1/15/19.
//

#include <stuka/dLP/naive_benders_subproblem.h>

stuka::dLP::NaiveBendersSubproblem::NaiveBendersSubproblem(const stuka::dLP::Subproblem sub, const Options opts) : sub_(sub), opts_(opts) {
  current_value_ = INF;

  prog_.c = sub.c;
  prog_.A_ub = sub.A_ub;
  prog_.b_ub = sub.b_ub;
  prog_.A_eq = sub.A_eq;
  prog_.b_eq = sub.b_eq;
  prog_.lb = sub.lb;
  prog_.ub = sub.ub;
}

stuka::dLP::NaiveBendersSubproblem::NaiveBendersSubproblem(stuka::dLP::NaiveBendersSubproblem &&sub) noexcept :
    NaiveBendersSubproblem(sub.sub_, sub.opts_) {}

stuka::dLP::BendersCut stuka::dLP::NaiveBendersSubproblem::getBendersCut(const Eigen::VectorXd &x) {

  if (sub_.b_ub) prog_.b_ub = std::make_shared<Eigen::VectorXd>(*sub_.b_ub - *sub_.C_ub * x);
  if (sub_.b_eq) prog_.b_eq = std::make_shared<Eigen::VectorXd>(*sub_.b_eq - *sub_.C_eq * x);

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
